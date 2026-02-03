use std::path::PathBuf;
use noodles_sam::Header;
use noodles_util::alignment as noodles_alignment;
use noodles_core::Region;

pub struct Alignment{
    reader: noodles_alignment::io::IndexedReader<std::fs::File>,
    header: noodles_sam::Header,
}

impl Alignment{

    pub fn new(alignment_path: PathBuf) -> Result<Self, Box<dyn std::error::Error>>{
        let mut builder = noodles_alignment::io::indexed_reader::Builder::default();
        let mut reader = builder.build_from_path(&alignment_path)?;
        let header=  reader.read_header()?;
        Ok(Alignment { reader, header })
    }

    pub fn get_region_coverage(&mut self, region: String) -> Result<i32, Box<dyn std::error::Error>>{
        let region_parsed: Region = region.parse()?;
        let query = self.reader.query(&self.header, &region_parsed)?;
        let mut count = 0;
        for _ in query{
            count += 1;
        }

        Ok(count)
    }
}