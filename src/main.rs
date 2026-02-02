mod utils;

use clap::Parser;
use noodles_sam::header::record::value::map::Inner;
use std::{env, path::{Path, PathBuf}};
use crate::utils::bam_handler::{get_chromosome_names, 
    get_chromosome_size_map, 
    bam_index_reader,
    get_chr_chunk_reads};
use crate::utils::alignment_handler;
use std::any::type_name;


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to the bam file
    #[arg(short, long)]
    bam_file_path: PathBuf,
    #[arg(long, default_value_t = 50)]
    bin_size: i32,
}


fn type_of<T>(_: T) -> &'static str {
    type_name::<T>()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    println!("{:?}", args);
    let bam_file_path = args.bam_file_path;
    let bam_index_file_name: PathBuf = PathBuf::from(format!("{}.bai", bam_file_path.display()));
    println!("{:#?}", bam_index_file_name);
    let bam_index = bam_index_reader(bam_index_file_name).unwrap();
    //println!("{:#?}", bam_index);
    let chunk: String = String::from("chr1:10000-20000");
    let bin_size: usize = args.bin_size as usize;
    let chromosome_size_map = get_chromosome_size_map(&bam_file_path).unwrap();
    let mut alignment = alignment_handler::Alignment::new(
        bam_file_path
    )?;
    
    for (chromosome, ref_seq) in &chromosome_size_map {
        let chr_len = ref_seq.length().get();
        let mut start = 1;

        while start < chr_len {
            let end = std::cmp::min(start + bin_size, chr_len);
            
            let region = format!("{}:{}-{}", chromosome, start, end);
            println!("{:#?}", region);
            println!("{:#?}", alignment.get_region_coverage(region)?);
            //println!("{:#?}", get_chr_chunk_reads(&bam_file_path, region).unwrap());
            start += bin_size;
        }
        break;
    }
    
    // println!("{:#?}", get_chromosome_names(args.bam_file_name).unwrap());
    // println!("{:#?}", get_chromosome_size_map(args.bam_file_name).unwrap());
    Ok(())
}