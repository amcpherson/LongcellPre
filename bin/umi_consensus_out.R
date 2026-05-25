#!/usr/bin/env Rscript

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Create argument parser
parser_args <- list()
for (i in seq(1, length(args), by = 2)) {
  if (i < length(args)) {
    param_name <- gsub("^--", "", args[i])
    param_value <- args[i + 1]
    parser_args[[param_name]] <- param_value
  }
}

# Extract required arguments
umi_count_path <- parser_args[["umi_count_path"]]
out_dir <- parser_args[["out_dir"]]
genome_path <- parser_args[["genome_path"]]
genome_name <- parser_args[["genome_name"]]
gene_bed_path <- parser_args[["gene_bed_path"]]

# Extract optional arguments with defaults
minimap2 <- parser_args[["minimap2"]] %||% "minimap2"
samtools <- parser_args[["samtools"]] %||% "samtools"
minimap_bed_path <- parser_args[["minimap_bed_path"]]
cores <- as.numeric(parser_args[["cores"]] %||% "1")
force_fastq_out <- parser_args[["force_fastq_out"]] %||% "FALSE"
force_fastq_out <- as.logical(force_fastq_out)
force_map <- parser_args[["force_map"]] %||% "FALSE"
force_map <- as.logical(force_map)

# Increase future globals size limit for parallel workers
options(future.globals.maxSize = 24 * 1024^3)

# Load libraries
library(LongcellPre)
library(Rsamtools)

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load input data
cat("Loading UMI count data...\n")
umi_count <- read.table(umi_count_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
gene_bed <- readRDS(gene_bed_path)

# Run UMI consensus output (fastq generation + remapping)
cat("Generating UMI-collapsed fastq and remapping...\n")
UMI_consensus_out(
  umi_count = umi_count,
  dir = out_dir,
  genome_path = genome_path,
  genome_name = genome_name,
  gene_bed = gene_bed,
  minimap2 = minimap2,
  samtools = samtools,
  minimap_bed_path = minimap_bed_path,
  cores = cores,
  force_fastq_out = force_fastq_out,
  force_map = force_map
)

cat("UMI consensus output completed successfully.\n")
