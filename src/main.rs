mod utils;

use clap::Parser;
use std::{env, path::{Path, PathBuf}};
use crate::utils::bam_handler::bam_index_reader;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to the bam file
    #[arg(short, long)]
    bam_file_name: PathBuf,
}

fn main() {
    let args = Cli::parse();
    println!("{:?}", args);
    let mut bam_index_file_name: PathBuf = PathBuf::from(format!("{}.bai", args.bam_file_name.display()));
    println!("{:#?}", bam_index_file_name);
    let bam_index = bam_index_reader(bam_index_file_name).unwrap();
    println!("{:#?}", bam_index);
}