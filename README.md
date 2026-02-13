# bamcowig

Converts BAM/CRAM alignment files to BigWig format with binned coverage.

I vibe code a lot, BUT REPO is written by me. Thats why it is going slow. 

## What it does

- Reads BAM or CRAM files and calculates read coverage across the genome in fixed-size bins
- Filters reads by mapping quality, alignment flags, strand, and duplicate status
- Supports normalization: CPM, RPKM, RPGC, BPM
- Handles both single-end and paired-end reads
- Outputs a BigWig file for genome browser visualization
- Processes chromosomes in parallel

## Build

```
cargo build --release
```

## Usage

```
bamcowig --bam-file-path <BAM_FILE/CRAM_FILE> --index-file-path <INDEX_FILE(bai/csi/crai)> --output-file <OUTPUT.bw>
```

### Options

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--bam-file-path` | `-b` | required | Path to BAM/CRAM file |
| `--index-file-path` | `-i` | required | Path to index file (.bai, .csi, .crai) |
| `--output-file` | `-o` | `coverage_over_bins.bed` | Output BigWig file |
| `--bin-size` | | `50` | Bin size in base pairs |
| `--threads` | `-t` | `8` | Number of threads |
| `--normalize` | | `none` | Normalization method: none, cpm, rpkm, rpgc, bpm |
| `--extend-to-fragment` | | `false` | Extend reads to fragment size using template length |
| `--fraction-counts` | `-f` | `false` | Pro-rate coverage for partial bin overlaps |

### Examples

```bash
# basic conversion
bamcowig -b sample.bam -i sample.bai -o sample.bw

# custom bin size and normalization
bamcowig -b sample.bam -i sample.bai -o sample.bw --bin-size 100 --normalize cpm

# paired-end with fragment extension
bamcowig -b sample.bam -i sample.bai -o sample.bw --extend-to-fragment --fraction-counts
```
