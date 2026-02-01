use std::{i32, path::PathBuf, path::Path, vec};
use noodles_sam::{self as sam, header::record::value::{Map, map::{ReferenceSequence, reference_sequence}}};
use noodles_bam as bam;

struct Alignment {
    alignment_file_path: PathBuf,
    header: sam::Header,
    // Fields computed later
    coverage: Option<Vec<f64>>,
    read_count: Option<i32>,
};

impl Alignment {
    fn new(alignment_file_path: PathBuf) -> Result<Alignment, Box<dyn std::error::Error>> {
        let mut reader = bam::io::reader::Builder::default().build_from_path(alignment_file_path.clone())?;
        let header = reader.read_header()?;
        Ok(Alignment {
            alignment_file_path,
            header,
            coverage: None,      // Not set yet
            read_count: None,     // Not set yet
        })
    }
}
