mod utils;

use clap::Parser;
use noodles_sam::header::record::value::map::Inner;
use crate::utils::filter::Filter;
use std::{env, path::{Path, PathBuf}};
use crate::utils::{bam_handler::{bam_index_reader, get_chr_chunk_reads, get_chromosome_names, get_chromosome_size_map}, filter};
use crate::utils::alignment_handler;
use std::any::type_name;
use std::fs;
 
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to the bam file
    #[arg(short, long)]
    bam_file_path: PathBuf,
    #[arg(short, long)]
    index_file_path: PathBuf,
    #[arg(long, default_value_t = 50)]
    bin_size: u16,
    #[arg(short, long, default_value = "coverage_over_bins.bed")]
    output_file: PathBuf,
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
}


fn type_of<T>(_: T) -> &'static str {
    type_name::<T>()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    println!("{:?}", args);
    let bam_file_path = args.bam_file_path;
    //let bam_index_file_name: PathBuf = PathBuf::from(format!("{}.bai", bam_file_path.display()));
    let bam_index_file = args.index_file_path;
    let bin_size = args.bin_size;
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let filter = Filter::default();
    let chromosome_size_map = get_chromosome_size_map(&bam_file_path).unwrap();
    let mut alignment = alignment_handler::Alignment::from_bam(
        bam_file_path,  bam_index_file
    )?;
    let coverage_over_bins_all_chromosomes = alignment.coverage_by_bin_all(bin_size, filter)?;
    write_output(args.output_file, coverage_over_bins_all_chromosomes, bin_size as usize)?;
    // println!("{:#?}", get_chromosome_names(args.bam_file_name).unwrap());
    // println!("{:#?}", get_chromosome_size_map(args.bam_file_name).unwrap());
    Ok(())
}

fn write_output(output: PathBuf, coverage_over_bins_all_chromosomes: Vec<Vec<f64>>, bin_size: usize) -> Result<(), Box<dyn std::error::Error>>{
    let mut content = String::new();
    for (chr_idx, coverage_over_bins) in coverage_over_bins_all_chromosomes.iter().enumerate() {
        for (bin_idx, &coverage_single_bin) in coverage_over_bins.iter().enumerate() {
            let bin_start = bin_idx * bin_size;
            let bin_end = bin_start + bin_size - 1;
            content.push_str(&format!("chr{}\t{}\t{}\n", chr_idx, bin_idx, coverage_single_bin));
        }
    }
    fs::write(output, content).expect("Should be able to write to `/foo/tmp`");
    Ok(())
}