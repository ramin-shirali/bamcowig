


use noodles_util::alignment::Record;
use noodles_sam::alignment::Record as _;

#[derive(Default)]
pub enum PairFilter {
    Strict,    // Only properly paired
    #[default]
    Lenient,   // Keep if read is mapped (even if mate unmapped)
    Off,       // No pair filtering at all
}


#[derive(PartialEq)]
pub enum StrandSelection {
    Both,
    Forward,
    Reverse,
}

pub struct Filter{
    minimum_mapping_quality: i32,
    ignore_duplicates_flag: bool,
    strand_selection: StrandSelection,
    secondary_alignment_skip: bool,
    supplementary_alignment_skip: bool,
    pair_filter: PairFilter,
}


impl Default for Filter{
    fn default() -> Filter {
        Filter {minimum_mapping_quality: 10, 
            ignore_duplicates_flag:false, 
            strand_selection: StrandSelection::Both,
            secondary_alignment_skip: true,
            supplementary_alignment_skip: true,
            pair_filter: PairFilter::Strict,
        }
   }
}

impl Filter{
    pub fn apply(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>> {
        if self.check_alignment(record)?{
            return Ok(true);
        }
        if self.minimum_mapping_quality > 0 && self.check_mapping_quality(record)? {
            return Ok(true);
        }
        if !self.ignore_duplicates_flag && self.check_duplicate(record)? {
            return Ok(true);
        }
        if self.strand_selection != StrandSelection::Both && self.check_strand_selection(record)? {
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
        if quality == 255{
            Ok(false) // 255 means quality score is not available
        }else{
            Ok(quality < self.minimum_mapping_quality)
        }
    }
    fn check_duplicate(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        let duplicate_alignment = record.flags()?.is_duplicate();
        Ok(duplicate_alignment)
    }
    fn check_strand_selection(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        let flags = record.flags()?;
        match self.strand_selection {
            StrandSelection::Both => Ok(false),
            StrandSelection::Forward => Ok(flags.is_reverse_complemented()), // Pay attention - true means the read is filtered. 
            StrandSelection::Reverse => Ok(!flags.is_reverse_complemented()), // Thats why it is reversed.
        }
    }
    fn check_secondary_alignment(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        let secondary_alignment = record.flags()?.is_secondary();
        Ok(secondary_alignment)
    }
    fn check_supplementary_alignment(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>>{
        let supplementary_alignment = record.flags()?.is_supplementary();
        Ok(supplementary_alignment)
    }
    fn check_alignment(&self, record: &Record) -> Result<bool, Box<dyn std::error::Error>> {
        let flags = record.flags()?;

        match self.pair_filter {
            PairFilter::Off => {
                Ok(false)
            }
            PairFilter::Lenient => {
                // Keep if read itself is mapped, ignore mate status
                Ok(flags.is_unmapped())
            }
            PairFilter::Strict => {
                if flags.is_segmented() {
                    // Paired-end: must be properly paired
                    Ok(!flags.is_properly_segmented())
                } else {
                    // Single-end: must be mapped
                    Ok(flags.is_unmapped())
                }
            }
        }
    }
}