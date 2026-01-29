use std::path::PathBuf;
use noodles_bam as bam;
use noodles_sam as sam;

pub fn bam_header_reader(input_path: PathBuf) -> Result<sam::Header, std::io::Error>{
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