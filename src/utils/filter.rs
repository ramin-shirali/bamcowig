

use std::io::Error;

use noodles_util::alignment::Record;
use noodles_sam::alignment::Record as _;

pub struct Filter{
    minimum_mapping_quality: i32,
    ignore_duplicates_flag: bool,
    strand_selection: bool,
    secondary_alignment_skip: bool,
    supplementary_alignment_skip: bool,
}


impl Default for Filter{
    fn default() -> Filter {
        Filter {minimum_mapping_quality: 10, 
            ignore_duplicates_flag:false, 
            strand_selection: true,
            secondary_alignment_skip: true,
            supplementary_alignment_skip: true
            }
   }
}

impl Filter{
    pub fn apply(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>> {
        if self.minimum_mapping_quality > 0 && self.check_mapping_quality(record)? {
            return Ok(true);
        }
        if !self.ignore_duplicates_flag && self.check_duplicate(record)? {
            return Ok(true);
        }
        if self.strand_selection && self.check_strand_selection(record)? {
            return Ok(true);
        }
        if self.secondary_alignment_skip && self.check_secondary_alignment(record)? {
            return Ok(true);
        }
        if self.supplementary_alignment_skip && self.check_supplementary_alignment(record)? {
            return Ok(true);
        }
        Ok(false)
    }

    fn check_mapping_quality(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        let quality = record.mapping_quality()
        .transpose()? //flips Option and Result from mapping_quality
        .map(|mq| mq.get() as i32) // MappingQuality is a wrapper around u8, .get() extracts the u8, then cast to i32
        .unwrap_or(0); // None means unmap. So 0 in that case to filter.         

        Ok(quality >= self.minimum_mapping_quality)
    }
    fn check_duplicate(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        Ok(true)
    }
    fn check_strand_selection(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        Ok(true)
    }
    fn check_secondary_alignment(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        Ok(true)
    }
    fn check_supplementary_alignment(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        Ok(true)
    }
}