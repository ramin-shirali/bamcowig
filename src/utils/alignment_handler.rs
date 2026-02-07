
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

pub type CsiIndex = csi::binning_index::Index<IndexMap<usize, VirtualPosition>>;

pub enum CountableIndex {
    Bai(bai::Index),
    Csi(CsiIndex),
}

pub trait AlignmentIndex: Sized{
    fn load(index_path: PathBuf) -> Result<Self, Box<dyn std::error::Error>>;
    fn count_total_reads(&self) -> Result<Option<u64>,Box<dyn std::error::Error>>;
}

impl AlignmentIndex for CountableIndex {
    fn load(index_path: PathBuf) -> Result<Self, Box<dyn std::error::Error>>{
        if index_path.exists(){          
            match index_path.extension().and_then(|e| e.to_str()) {
                Some("bai") => {
                    let index = bai::fs::read(&index_path)?; 
                    return Ok(CountableIndex::Bai(index))
                }
                Some("csi") => {
                    let index = csi::fs::read(&index_path)?; 
                    return Ok(CountableIndex::Csi(index))
                }
                _ => return Err("Unsupported index format".into()),
            }
        }
        Err("No index (.bai or .csi) found".into())
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
    fn load(index_path: PathBuf) -> Result<Self, Box<dyn std::error::Error>>{
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




pub struct Alignment<I>{
    reader: noodles_alignment::io::IndexedReader<std::fs::File>,
    header: noodles_sam::Header,
    index: I,
    total_reads: u64,
    file_type: String
}


impl Alignment<CountableIndex> {
    pub fn from_bam(alignment_path: PathBuf, index_path: PathBuf) -> Result<Self, Box<dyn std::error::Error>> {
        let mut reader = noodles_alignment::io::indexed_reader::Builder::default()
            .build_from_path(&alignment_path)?;
        let header = reader.read_header()?;
        let index = CountableIndex::load(index_path)?;
        let total_reads = index.count_total_reads()?.unwrap_or(0);
        Ok(Alignment {
            reader,
            header,
            index,
            total_reads,
            file_type: "bam".to_string(),
        })
        
    }
}

impl Alignment<crai::Index> {
    pub fn from_cram(alignment_path: PathBuf, index_path: PathBuf) -> Result<Self, Box<dyn std::error::Error>> {
        let mut reader = noodles_alignment::io::indexed_reader::Builder::default()
            .build_from_path(&alignment_path)?;
        let header = reader.read_header()?;
        let index = crai::fs::read(index_path)?;
        let total_reads: u64 = Self::count_from_containers(&alignment_path)?;
        
        Ok(Alignment {
            reader,
            header,
            index,
            total_reads,
            file_type: "cram".to_string(),
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

    pub fn coverage_by_bin_all(&mut self, bin_size:u16) -> Result<Vec<Vec<u32>>, Box< dyn std::error::Error>>{
        let bin_size = bin_size as usize;
        let mut coverage_over_bins_per_chromosome: Vec<Vec<u32>> = Vec::new();
        let refs: Vec<_> = self.header.reference_sequences()
            .iter()
            .map(|(chr, info)| (chr.clone(), info.length().get()))
            .collect();

        for (chromosome, chromosome_length) in refs{
            let mut coverage_over_bins = self.get_coverage_chr(bin_size, chromosome.to_string(), chromosome_length)?;
            coverage_over_bins_per_chromosome.push(coverage_over_bins);
        }
        Ok(coverage_over_bins_per_chromosome)
    }

    pub fn get_coverage_chr(&mut self, bin_size: usize, chromosome: String, chromosome_length: usize) -> Result<Vec<u32>, Box<dyn std::error::Error>>{   
        let mut coverage_over_bins: Vec<u32> = Vec::new();
        let mut start = 0;
        while start < chromosome_length {
            let end = std::cmp::min(start + bin_size, chromosome_length);
            
            let region = format!("{}:{}-{}", chromosome, start, end);
            coverage_over_bins.push(self.get_region_coverage(region)?);
            //println!("{:#?}", get_chr_chunk_reads(&bam_file_path, region).unwrap());
            start += bin_size;
        }
        Ok(coverage_over_bins)
    }

    pub fn get_region_coverage(&mut self, region: String) -> Result<u32, Box<dyn std::error::Error>>{
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