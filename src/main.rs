mod utils;

use clap::Parser;
use std::{env, path::{Path, PathBuf}};
use crate::utils::bam_handler::{get_chromosome_names, 
    get_chromosome_size_map, 
    bam_index_reader,
    get_chr_chunk_reads};
use std::any::type_name;


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to the bam file
    #[arg(short, long)]
    bam_file_name: PathBuf,
}


fn type_of<T>(_: T) -> &'static str {
    type_name::<T>()
}

fn main() {
    let args = Cli::parse();
    println!("{:?}", args);
    let bam_index_file_name: PathBuf = PathBuf::from(format!("{}.bai", args.bam_file_name.display()));
    println!("{:#?}", bam_index_file_name);
    let bam_index = bam_index_reader(bam_index_file_name).unwrap();
    //println!("{:#?}", bam_index);
    let chunk: String = String::from("chr1:10000-20000");
    println!("{:#?}", get_chr_chunk_reads(args.bam_file_name, chunk));
    // println!("{:#?}", get_chromosome_names(args.bam_file_name).unwrap());
    // println!("{:#?}", get_chromosome_size_map(args.bam_file_name).unwrap());
}