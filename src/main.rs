use clap::Parser;
use std::{env, path::PathBuf};

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
}