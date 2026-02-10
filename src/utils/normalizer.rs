
use rayon::prelude::*;

pub fn cpm(coverage_over_bins_all_chromosomes: Vec<Vec<f64>>, total_read_count: u64) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error + Send + Sync>>{
    
    Ok(
        coverage_over_bins_all_chromosomes.par_iter()
        .map(|chr| {
            chr.iter()
                .map(|&count| count as f64 * 1_000_000.0 / total_read_count as f64)
                .collect()
        })
        .collect()
    )
    
}


pub fn rpkm(coverage_over_bins_all_chromosomes: Vec<Vec<f64>>, total_read_count: u64, bin_size: usize) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error + Send + Sync>>{
    Ok(
        coverage_over_bins_all_chromosomes.par_iter()
        .map(|chr| {
            chr.iter()
                .map(|&count| (count as f64 * 1_000_000_000.0) / (total_read_count as f64 * bin_size as f64))
                .collect()
        })
        .collect()
    )
}


pub fn rpgc(coverage_over_bins_all_chromosomes: Vec<Vec<f64>>, total_read_count: u64, effective_genome_size: usize, average_read_length: usize) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error + Send + Sync>>{
    let scale = effective_genome_size as f64 / (total_read_count as f64 * average_read_length as f64);
    Ok(
        coverage_over_bins_all_chromosomes.par_iter()
        .map(|chr| {
            chr.iter()
                .map(|&count| count as f64 * scale)
                .collect()
        })
        .collect()
    )
}


pub fn bpm(coverage_over_bins_all_chromosomes: Vec<Vec<f64>>) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error + Send + Sync>>{
    let total_bins_count = total_bins_count(&coverage_over_bins_all_chromosomes)?;
    Ok(    
        coverage_over_bins_all_chromosomes.par_iter()
        .map(|chr| {
            chr.iter()
                .map(|&count| count as f64 * 1_000_000.0 / total_bins_count as f64)
                .collect()
        })
        .collect()
    )
}

fn total_bins_count(coverage_over_bins_all_chromosomes: &Vec<Vec<f64>>) -> Result<f64, Box<dyn std::error::Error + Send + Sync>>{

    Ok(
        coverage_over_bins_all_chromosomes
        .iter()
        .flatten()
        .sum()
    )

}
