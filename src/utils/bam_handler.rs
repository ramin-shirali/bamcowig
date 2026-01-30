use std::{path::PathBuf, vec};
use noodles_bam as bam;
use noodles_sam as sam;
use bstr;

fn bam_header_reader(input_path: PathBuf) -> Result<sam::Header, std::io::Error>{
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

pub fn get_chromosome_names(input_path: PathBuf) -> Result<Vec<bstr::BString>, std::io::Error>{
    let header: sam::Header = bam_header_reader(input_path)?;
    let chromosomes = header.reference_sequences().keys().cloned().collect();
    // Ok::<_, io::Error>(())
    Ok(chromosomes)
}