use std::path::PathBuf;

use noodles_cram as cram;
use noodles_sam as sam;

pub fn cram_header_reader(input_path: PathBuf) -> Result<sam::Header, std::io::Error>{
    let mut reader = cram::io::reader::Builder::default().build_from_path(input_path)?;
    let header = reader.read_header()?;
    Ok(header)
    // Ok::<_, io::Error>(())
}