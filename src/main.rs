mod utils;

use clap::Parser;
use noodles_sam::header::record::value::map::Inner;
use crate::utils::filter::Filter;
use std::{collections::HashMap, env, path::{Path, PathBuf}};
use crate::utils::alignment_handler;
use std::any::type_name;
use std::fs;
use crate::utils::normalizer::{cpm, rpkm, rpgc, bpm};
use bigtools::{BigWigWrite, Value};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::bed::bedparser::BedIteratorStream;

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
    #[arg(long, default_value_t = false)]
    extend_to_fragment: bool,
    #[arg(long, default_value = "none")]
    normalize: String,
    #[arg(short, long, default_value_t = false)]
    fraction_counts: bool,
}


fn type_of<T>(_: T) -> &'static str {
    type_name::<T>()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    println!("{:?}", args);
    let bam_file_path = args.bam_file_path;
    
    let bam_index_file = args.index_file_path;
    let bin_size = args.bin_size;
    let extend_to_fragment = args.extend_to_fragment;
    let fraction_counts = args.fraction_counts;
    let max_threads = args.threads;

    rayon::ThreadPoolBuilder::new()
        .num_threads(max_threads)
        .build_global()
        .unwrap();

    let filter = Filter::default();
    let mut alignment = alignment_handler::Alignment::from_bam(
        bam_file_path,  bam_index_file, None
    )?;

    let chromosome_names = alignment.get_chromosome_names_str()?;
    let chromosome_sizes = alignment.get_chromosome_sizes()?;

    let coverage_over_bins_all_chromosomes = alignment.coverage_by_bin_all(bin_size, filter, extend_to_fragment, fraction_counts)?;
    let normalized_over_bins_all_chromosomes = cpm(coverage_over_bins_all_chromosomes, alignment.total_reads().clone()).unwrap();
    write_bigwig_output(args.output_file, normalized_over_bins_all_chromosomes, chromosome_names, chromosome_sizes, bin_size as usize, max_threads)?;
    
    Ok(())
}


fn write_bigwig_output(output: PathBuf, coverage_over_bins_all_chromosomes: Vec<Vec<f64>>, chrom_names: Vec<String>, chrom_sizes: Vec<usize>, bin_size: usize, threads: usize) -> Result<(), Box<dyn std::error::Error>>{
    
    let runtime = tokio::runtime::Builder::new_multi_thread()
        .worker_threads(threads)
        .enable_all()
        .build()?;

    let mut sorted: Vec<_> = chrom_names
        .into_iter()
        .zip(chrom_sizes)
        .zip(coverage_over_bins_all_chromosomes)
        .map(|((name, size), coverage)| (name, size, coverage))
        .collect();
    
    sorted.sort_by(|a, b| a.0.cmp(&b.0));

    let chrom_map: HashMap<String, u32> = sorted
        .iter()
        .map(|(name, size, _)| (name.to_string(), *size as u32))
        .collect();

    let values_iter = sorted
        .into_iter()
        .flat_map(|(chrom_name, _, bins)|
        {
            bins.into_iter()
                .enumerate()
                .filter(|(_, val)| *val != 0.0)
                .map(move |(bin_idx, val)| {
                    let start = (bin_idx as u32) * bin_size as u32;
                    let end = start + bin_size as u32;
                    (chrom_name.clone(), Value { start, end, value: val as f32 })
                })
        });
    

    let data_source = BedParserStreamingIterator::wrap_infallible_iter(values_iter, false);
    let writer = BigWigWrite::create_file(output.to_string_lossy().to_string(), chrom_map)?;
    writer.write(data_source, runtime)?;

    Ok(())
}

fn write_bed_output(output: PathBuf, coverage_over_bins_all_chromosomes: Vec<Vec<f64>>, bin_size: usize) -> Result<(), Box<dyn std::error::Error>>{
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