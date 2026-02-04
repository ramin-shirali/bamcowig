

use std::io::Error;

use noodles_util::alignment::Record;

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
            ignore_duplicates_flag:true, 
            strand_selection: true,
            secondary_alignment_skip: true,
            supplementary_alignment_skip: true
            }
   }
}

impl Filter{
    pub fn apply(record: Record) -> Result<Record, Box<dyn std::error::Error>>{
        Ok(record)
    }
}