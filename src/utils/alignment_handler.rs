use std::{i32, path::PathBuf, path::Path, vec, fs::File, ffi::OsStr};
use noodles_sam::{self as sam, header::record::value::{Map, map::{ReferenceSequence, reference_sequence}}};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_bgzf;
use noodles_core::Region;

type BamIndexedReader = bam::io::IndexedReader<noodles_bgzf::io::reader::Reader<File>>;
type CramIndexedReader = cram::io::IndexedReader<File>;


enum AlignmentReader {
    Bam(BamIndexedReader),
    Cram(CramIndexedReader),
}


pub struct Alignment {
    alignment_file_path: PathBuf,
    reader: AlignmentReader,
    header: sam::Header,
    // Fields computed later
    coverage: Option<Vec<f64>>,
    read_count: Option<i32>,
}

impl Alignment {
    pub fn new(alignment_file_path: PathBuf) -> Result<Alignment, Box<dyn std::error::Error>> {
        let input_alignment_extension = alignment_file_path.extension().and_then(OsStr::to_str).unwrap_or("");
        let (reader, header) = match input_alignment_extension {
            "bam" => {
                let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(alignment_file_path.clone())?;
                let header = reader.read_header()?;
                (AlignmentReader::Bam(reader), header)
            }
            "cram" => {
                let mut reader = cram::io::indexed_reader::Builder::default().build_from_path(alignment_file_path.clone())?;
                let header = reader.read_header()?;
                (AlignmentReader::Cram(reader), header)
            }
            _ => return Err("Unsupported format. Use .bam or .cram".into()),
        };

        Ok(Alignment {
            alignment_file_path,
            reader,
            header,
            coverage: None,      // Not set yet
            read_count: None,     // Not set yet
        })
    }

    

    pub fn get_region_coverage(&mut self, region: String) -> Result<i32, Box<dyn std::error::Error>>{

        let header = &self.header;
        let mut count = match &mut self.reader{
            AlignmentReader::Bam(reader) => {
                Self::get_region_coverage_bam(header, reader, region)
            }
            AlignmentReader::Cram(reader) => {
                Self::get_region_coverage_cram(header, reader, region)
            }
        };
    
        Ok(count.unwrap())
    }

    fn get_region_coverage_bam(header: &sam::Header, reader: &mut BamIndexedReader, region:String) -> Result<i32, Box<dyn std::error::Error>>{
        let region = region.parse().unwrap();
        let query = reader.query(header, &region)?;
        let mut count = 0;
        query.records().inspect(|_| count += 1).collect::<Vec<_>>();
        Ok(count)
    }

    fn get_region_coverage_cram(header: &sam::Header, reader: &mut CramIndexedReader, region:String, ) -> Result<i32, Box<dyn std::error::Error>>{
        let region = region.parse().unwrap();
        let query = reader.query(header, &region)?;
        let mut count = 0;
        query.inspect(|_| count += 1).collect::<Vec<_>>();
        Ok(count)
    }
    
}
