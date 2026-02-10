
use std::path::PathBuf;
use noodles_bam::bai;
use noodles_bgzf::VirtualPosition;
use noodles_util::alignment as noodles_alignment;
use noodles_core::Region;
use noodles_cram::{crai, io::reader};
use noodles_csi as csi;
use noodles_csi::BinningIndex;
use indexmap::IndexMap;
use csi::binning_index::index::reference_sequence::Index as ReferenceSequenceIndex;
use noodles_cram::io::reader as CramReader;
use noodles_cram::io::reader::Container as CramContainer;
use rayon::{prelude::*};
use std::path::Path;
use crate::filter::Filter;
pub type CsiIndex = csi::binning_index::Index<IndexMap<usize, VirtualPosition>>;
use getset::{Getters, Setters, MutGetters};

pub enum CountableIndex {
    Bai(bai::Index),
    Csi(CsiIndex),
}

pub trait AlignmentIndex: Sized{
    fn load(index_path: &Path) -> Result<Self, Box<dyn std::error::Error>>;
    fn count_total_reads(&self) -> Result<Option<u64>,Box<dyn std::error::Error>>;
}

impl AlignmentIndex for CountableIndex {


    fn load(index_path: &Path) -> Result<Self, Box<dyn std::error::Error>>{
        if !index_path.exists() {
            return Err("No index (.bai or .csi) found".into());
        }         
        match index_path.extension().and_then(|e| e.to_str()) {
            Some("bai") => {
                let index = bai::fs::read(index_path)?; 
                Ok(CountableIndex::Bai(index))
            }
            Some("csi") => {
                let index = csi::fs::read(index_path)?; 
                Ok(CountableIndex::Csi(index))
            }
            _ => Err("Unsupported index format".into()),
        }

    }
    
    fn count_total_reads(&self) -> Result<Option<u64>, Box<dyn std::error::Error>> {
        match self {
            CountableIndex::Bai(index) => count_from_index(index),
            CountableIndex::Csi(index) => count_from_index(index),
        }
    }
}

fn count_from_index<I>(index: &csi::binning_index::Index<I>) -> 
    Result<Option<u64>, Box<dyn std::error::Error>>
    where I: ReferenceSequenceIndex
    {
        use csi::binning_index::index::reference_sequence::Bin;

        let depth = index.depth();
        let metadata_bin_id = Bin::metadata_id(depth);
        let mut count: u64 = 0;
        
        for reference in index.reference_sequences() {
            let bins = reference.bins(); // Look for the metadata bin 
            if let Some(metadata_bin) = bins.get(&metadata_bin_id) { 
                let chunks = metadata_bin.chunks(); // Metadata bin has exactly 2 chunks: // chunks[0].end = mapped count (stored as virtual position) // chunks[1].end = unmapped count 
                if chunks.len() >= 2 { // Virtual position encodes count in compressed() part 
                    let mapped = chunks[0].end().compressed(); 
                    let unmapped = chunks[1].end().compressed(); 
                    count += mapped; 
                    count += unmapped; 
                } 
            } 
        }
        if let Some(unmapped) = index.unplaced_unmapped_record_count() {
            count += unmapped;
        }
        Ok(Some(count))
    }


impl AlignmentIndex for crai::Index{
    fn load(index_path: &Path) -> Result<Self, Box<dyn std::error::Error>>{
        if index_path.exists(){
            let index = crai::fs::read(index_path)?;
            return Ok(index);
        }
        Err("No .crai index found".into())
    }
    fn count_total_reads(&self) -> Result<Option<u64>, Box<dyn std::error::Error>>{
        println!("count of reads is not in cram.crai file. use alignment.count_total_reads()");
        Ok(None)
    }    
}



#[derive(Getters, Setters, MutGetters)]
#[getset(get = "pub")]
pub struct Alignment<I>{
    file_path: PathBuf,
    index_path: PathBuf,
    reader: noodles_alignment::io::IndexedReader<std::fs::File>,
    header: noodles_sam::Header,
    index: I,
    total_reads: u64,
    file_type: String,
    is_pair_end: bool,
}


impl Alignment<CountableIndex> {
    pub fn from_bam(alignment_path: PathBuf, index_path: PathBuf, pair_end_flag: Option<bool>) -> Result<Self, Box<dyn std::error::Error>> {
        let mut reader = noodles_alignment::io::indexed_reader::Builder::default()
            .build_from_path(&alignment_path)?;
        let header = reader.read_header()?;
        let index = CountableIndex::load(&index_path)?;
        let total_reads = index.count_total_reads()?.unwrap_or(0);
        let is_pair_end = pair_end_flag.unwrap_or(reader.records(&header).next().unwrap()?.flags()?.is_segmented());
        Ok(Alignment {
            file_path: alignment_path,
            index_path: index_path,
            reader,
            header,
            index,
            total_reads,
            file_type: "bam".to_string(),
            is_pair_end,
        })
        
    }
}

impl Alignment<crai::Index> {
    pub fn from_cram(alignment_path: PathBuf, index_path: PathBuf, pair_end_flag: Option<bool>) -> Result<Self, Box<dyn std::error::Error>> {
        let mut reader = noodles_alignment::io::indexed_reader::Builder::default()
            .build_from_path(&alignment_path)?;
        let header = reader.read_header()?;
        let index = crai::fs::read(&index_path)?;
        let total_reads: u64 = Self::count_from_containers(&alignment_path)?;
        let is_pair_end = pair_end_flag.unwrap_or(reader.records(&header).next().unwrap()?.flags()?.is_segmented());
        Ok(Alignment {
            file_path: alignment_path,
            index_path,
            reader,
            header,
            index,
            total_reads,
            file_type: "cram".to_string(),
            is_pair_end,
        })
    }
}


impl<I> Alignment<I>{

    fn count_by_iteration(reader: &mut noodles_alignment::io::IndexedReader<std::fs::File>, header: &noodles_sam::Header) -> Result<u64, Box<dyn std::error::Error>>{
        let mut count  = 0u64;
        for result in reader.records(header){
            let _ = result?;
            count += 1;
        }   
        Ok(count)
    }

    pub fn coverage_by_bin_all(& mut self, bin_size:u16, filter: Filter, extend_to_fragment: bool) -> Result<Vec<Vec<f64>>, Box< dyn std::error::Error>>{
        let bin_size = bin_size as usize;

        let refs: Vec<_> = self.header.reference_sequences()
            .iter()
            .map(|(chr, info)| (chr.clone(), info.length().get()))
            .collect();

        let file_path = &self.file_path;
        if extend_to_fragment{ // /extend by fragments
            let coverage_over_bins_per_chromosome: Vec<Vec<f64>> = refs.par_iter()
            .map(|(chromosome, chromosome_length)|{
                let mut reader = noodles_alignment::io::indexed_reader::Builder::default()
                .build_from_path(file_path)?;
            reader.read_header()?;
            Self::get_coverage_chr_with_reader_iterating_reads_extend_to_fragment(
                &mut reader, &self.header, bin_size,
                chromosome.to_string(), *chromosome_length, filter.clone(), self.is_pair_end
            )
            }).collect::<Result<Vec<Vec<f64>>, _>>()
            .map_err(|e| e as Box<dyn std::error::Error>)?;
        
            Ok(coverage_over_bins_per_chromosome)
            
        }else{ // just calculate aligned regions for coverage
            let coverage_over_bins_per_chromosome: Vec<Vec<f64>> = refs.par_iter()
            .map(|(chromosome, chromosome_length)|{
                let mut reader = noodles_alignment::io::indexed_reader::Builder::default()
                .build_from_path(file_path)?;
            reader.read_header()?;
            Self::get_coverage_chr_with_reader_iterating_reads(
                &mut reader, &self.header, bin_size,
                chromosome.to_string(), *chromosome_length, filter.clone()
            )
            }).collect::<Result<Vec<Vec<f64>>, _>>()
            .map_err(|e| e as Box<dyn std::error::Error>)?;
        
            Ok(coverage_over_bins_per_chromosome)
        }
        
    }



    fn get_coverage_chr_with_reader_iterating_reads_extend_to_fragment(
        reader: &mut noodles_alignment::io::IndexedReader<std::fs::File>,
            header: &noodles_sam::Header,
            bin_size: usize,
            chromosome: String,
            chromosome_length: usize,
            filter: Filter,
            is_pair_end: bool,
        ) -> Result<Vec<f64>, Box<dyn std::error::Error + Send + Sync>> 
    {
        let bin_count = (chromosome_length / bin_size) +1 ;


        let region: Region = format!("{}:{}-{}", chromosome, 1, chromosome_length).parse()?; //single_chromosome
        if is_pair_end{
            Ok(Self::coverage_extend_to_fragment_pair_end(
            reader, 
            header, 
            bin_size, 
            filter,
            region,
            bin_count,)?)
        }else{
            Ok(Self::coverage_extend_to_fragment_single_end(
            reader, 
            header, 
            bin_size, 
            filter,
            region,
            bin_count,)?)
        }
        
    }

    fn coverage_extend_to_fragment_pair_end(
        reader: &mut noodles_alignment::io::IndexedReader<std::fs::File>,
            header: &noodles_sam::Header,
            bin_size: usize,
            filter: Filter,
            region: Region, 
            bin_count: usize,
        ) -> Result<Vec<f64>, Box< dyn std::error::Error + Send + Sync>>
    {
        let mut coverage_over_bins:Vec<f64> = vec![0f64; bin_count];
        reader.query(header, &region)?.map(|r| r.unwrap())
        .filter(|record| !filter.apply(record).unwrap_or(false)).filter(|record| record.template_length().unwrap() > 0)
        .for_each(|record| {

            let start = record.alignment_start().unwrap().unwrap();
            let template_length = record.template_length().unwrap();
            if template_length < 0{
                println!("sanity check!!! template length <0!! template_length: {:#?}", template_length)
            }
            let fragment_end = start.get() + template_length.abs() as usize;
            let start_bin = (start.get() - 1) /  bin_size; //noodles positions are 1-based. Yikes.
            let start_offset= (start.get() - 1) % bin_size;
            let end_bin = std::cmp::min((fragment_end - 1) / bin_size, bin_count - 1); // Some aligners (e.g. BWA) can produce alignments that extend past the reference end. The BAM spec doesn't enforce that. Yikes.

            let end_offset =  (fragment_end - 1) % bin_size;
            let flags = record.flags().unwrap();
            
            if end_bin == start_bin{
                coverage_over_bins[start_bin] += 1.0;
            }else{
                
                coverage_over_bins[start_bin] = (bin_size - start_offset) as f64 / bin_size as f64;
                coverage_over_bins[end_bin] = end_offset as f64 / bin_size as f64;
                if end_bin - start_bin > 1{
                    for i in (start_bin + 1)..end_bin{
                        coverage_over_bins[i] += 1.0;
                }
            }
            
            }
        });
        Ok(coverage_over_bins)
    }

    fn coverage_extend_to_fragment_single_end(
        reader: &mut noodles_alignment::io::IndexedReader<std::fs::File>,
            header: &noodles_sam::Header,
            bin_size: usize,
            filter: Filter,
            region: Region, 
            bin_count: usize,
        ) -> Result<Vec<f64>, Box< dyn std::error::Error + Send + Sync>>
    {
        let mut coverage_over_bins:Vec<f64> = vec![0f64; bin_count];
        reader.query(header, &region)?.map(|r| r.unwrap())
        .filter(|record| !filter.apply(record).unwrap_or(false))
        .for_each(|record| {

            let start = record.alignment_start().unwrap().unwrap().get();
            let end = record.alignment_end().unwrap().unwrap().get();
            let template_length = record.template_length().unwrap();
            let mut fragment_start: usize;
            let mut fragment_end: usize;

            if template_length < 0{
                fragment_start = (end as i32 + template_length) as usize;
                fragment_end = end;
            }else{
                fragment_start = start;
                fragment_end = start + template_length as usize;
            }
            
            let start_bin = (fragment_start - 1) /  bin_size; //noodles positions are 1-based. Yikes.
            let start_offset= (fragment_start - 1) % bin_size;
            let end_bin = std::cmp::min((fragment_end - 1) / bin_size, bin_count - 1); // Some aligners (e.g. BWA) can produce alignments that extend past the reference end. The BAM spec doesn't enforce that. Yikes.

            let end_offset =  (fragment_end - 1) % bin_size;
            
            if end_bin == start_bin{
                coverage_over_bins[start_bin] += 1.0;
            }else{
                
                coverage_over_bins[start_bin] = (bin_size - start_offset) as f64 / bin_size as f64;
                coverage_over_bins[end_bin] = end_offset as f64 / bin_size as f64;
                if end_bin - start_bin > 1{
                    for i in (start_bin + 1)..end_bin{
                        coverage_over_bins[i] += 1.0;
                }
            }
            
            }
        });
        Ok(coverage_over_bins)

        
    }

    fn get_coverage_chr_with_reader_iterating_reads(
        reader: &mut noodles_alignment::io::IndexedReader<std::fs::File>,
            header: &noodles_sam::Header,
            bin_size: usize,
            chromosome: String,
            chromosome_length: usize,
            filter: Filter,
        ) -> Result<Vec<f64>, Box<dyn std::error::Error + Send + Sync>> 
    {
        let bin_count = (chromosome_length / bin_size) +1 ;
        let mut coverage_over_bins:Vec<f64> = vec![0f64; bin_count];

        let region: Region = format!("{}:{}-{}", chromosome, 1, chromosome_length).parse()?; //single_chromosome

        reader.query(header, &region)?.map(|r| r.unwrap())
        .filter(|record| !filter.apply(record).unwrap_or(false))
        .for_each(|record| {
            let start = record.alignment_start().unwrap().unwrap();
            let end = record.alignment_end().unwrap().unwrap();
            let start_bin = (start.get() - 1) /  bin_size; //noodles positions are 1-based. Yikes.
            let start_offset= (start.get() - 1) % bin_size;
            let end_bin = std::cmp::min((end.get() - 1) / bin_size, bin_count - 1); // Some aligners (e.g. BWA) can produce alignments that extend past the reference end. The BAM spec doesn't enforce that. Yikes.

            let end_offset =  (end.get() - 1) % bin_size;
            let flags = record.flags().unwrap();
            
            if end_bin == start_bin{
                coverage_over_bins[start_bin] += 1.0;
            }else{
                //println!("bin_size {:#?} and offset {:#?} and start_bin {:#?} and bin_count {:#?}", bin_size, start_offset, start_bin, bin_count);

                coverage_over_bins[start_bin] = (bin_size - start_offset) as f64 / bin_size as f64;
                coverage_over_bins[end_bin] = end_offset as f64 / bin_size as f64;
                if end_bin - start_bin > 1{
                    for i in (start_bin + 1)..end_bin{
                        coverage_over_bins[i] += 1.0;
                }
            }
            
            }
        });
        Ok(coverage_over_bins)
    }


    fn get_coverage_chr_with_reader(
            reader: &mut noodles_alignment::io::IndexedReader<std::fs::File>,
            header: &noodles_sam::Header,
            bin_size: usize,
            chromosome: String,
            chromosome_length: usize,
        ) -> Result<Vec<u32>, Box<dyn std::error::Error + Send + Sync>> {
        let mut coverage_over_bins: Vec<u32> = Vec::new();
        let mut start = 1;
        while start < chromosome_length {
            let end = std::cmp::min(start + bin_size, chromosome_length);
            let region: Region = format!("{}:{}-{}", chromosome, start, end).parse()?;
            let query = reader.query(header, &region)?;
            let mut count = 0u32;
            for result in query {
                let _ = result?;
                count += 1;
            }
            coverage_over_bins.push(count);
            start += bin_size;
        }
        Ok(coverage_over_bins)
    }
    pub fn get_coverage_chr(&mut self, bin_size: usize, chromosome: String, chromosome_length: usize) -> Result<Vec<u32>, Box<dyn std::error::Error + Send + Sync>>{   
        let mut coverage_over_bins: Vec<u32> = Vec::new();
        let mut start = 1; //start position is 1. 0 gives error.
        while start < chromosome_length {
            let end = std::cmp::min(start + bin_size, chromosome_length);
            
            let region = format!("{}:{}-{}", chromosome, start, end);
            coverage_over_bins.push(self.get_region_coverage(region)?);
            //println!("{:#?}", get_chr_chunk_reads(&bam_file_path, region).unwrap());
            start += bin_size;
        }
        Ok(coverage_over_bins)
    }

    pub fn get_region_coverage(&mut self, region: String) -> Result<u32, Box<dyn std::error::Error + Send + Sync>>{
        let region_parsed: Region = region.parse()?;
        let query = self.reader.query(&self.header, &region_parsed)?;
        let mut count = 0u32;
        for result in query{
            let _ = result?;
            count += 1;
        }
        Ok(count)
        
    }

    fn count_from_containers(alignment_path: &PathBuf) -> Result<u64, Box<dyn std::error::Error>> {
    
        let mut reader = CramReader::Builder::default()
            .build_from_path(&alignment_path)?;

        // Required for CRAM
        reader.read_file_definition()?;
        reader.read_file_header()?;

        let mut total = 0u64;
        let mut container = CramContainer::default();

        while reader.read_container(&mut container)? > 0 {
            total += container.header().record_count() as u64;
        }

        Ok(total)
    }
}