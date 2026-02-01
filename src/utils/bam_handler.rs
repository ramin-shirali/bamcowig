use std::{i32, path::PathBuf, path::Path, vec};
use indexmap::IndexMap;
use noodles_bam as bam;
use noodles_sam::{self as sam, header::record::value::{Map, map::{ReferenceSequence, reference_sequence}}};
use bstr;

fn bam_header_reader(input_path: &Path) -> Result<sam::Header, std::io::Error>{
    let mut reader = bam::io::reader::Builder::default().build_from_path(input_path)?;
    let header = reader.read_header()?;
    Ok(header)
    // Ok::<_, io::Error>(())
}

pub fn bam_index_reader(input_path: PathBuf) -> Result<bam::bai::Index, std::io::Error>{
    let index = bam::bai::fs::read(input_path)?;
    // Ok::<_, io::Error>(())
    Ok(index)
}

pub fn get_chromosome_names(input_path: &Path) -> Result<Vec<bstr::BString>, std::io::Error>{
    let header: sam::Header = bam_header_reader(input_path)?;
    let chromosomes = header.reference_sequences().keys().cloned().collect();
    // Ok::<_, io::Error>(())
    Ok(chromosomes)
}

pub fn get_chromosome_size_map(input_path: &Path) -> Result<IndexMap<bstr::BString,Map<ReferenceSequence>>, std::io::Error>{
    let header: sam::Header = bam_header_reader(input_path)?;
    let reference_sequence = header.reference_sequences().clone();
    // Ok::<_, io::Error>(())
    Ok(reference_sequence)
}

// pub fn get_total_coverage(input_path: PathBuf) -> Result<IndexMap<bstr::BString,Map<ReferenceSequence>>, std::io::Error>{
//     let index = bam_index_reader(input_path).unwrap();
    
// }

pub fn get_chr_chunk_reads(input_path: &Path, chunk: String) -> Result<i32, std::io::Error>{
    let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(input_path)?;
    let header = reader.read_header()?;

    let region = chunk.parse().unwrap();
    let query = reader.query(&header, &region)?;

    // for result in query.records() {
    //     let record = result?;
    //     // ...
    // }
    let mut count = 0;
    let region_reads_vec = query.records().inspect(|_| count += 1).collect::<Vec<_>>();
    Ok(count)
    //Ok::<_, Box<dyn std::error::Error>>(())
}