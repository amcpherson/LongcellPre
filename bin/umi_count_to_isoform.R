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
gene_bed_path <- parser_args[["gene_bed_path"]]
gtf_path <- parser_args[["gtf_path"]]

# Extract optional arguments with defaults
gene_col <- parser_args[["gene_col"]] %||% "gene"
bed_gene_col <- parser_args[["bed_gene_col"]] %||% "gene"
bed_strand_col <- parser_args[["bed_strand_col"]] %||% "strand"
filter_only_intron <- parser_args[["filter_only_intron"]] %||% "TRUE"
filter_only_intron <- as.logical(filter_only_intron)
mid_offset_thresh <- as.numeric(parser_args[["mid_offset_thresh"]] %||% "3")
overlap_thresh <- as.numeric(parser_args[["overlap_thresh"]] %||% "0")
gtf_gene_col <- parser_args[["gtf_gene_col"]] %||% "gene"
gtf_start_col <- parser_args[["gtf_start_col"]] %||% "start"
gtf_end_col <- parser_args[["gtf_end_col"]] %||% "end"
gtf_iso_col <- parser_args[["gtf_iso_col"]] %||% "transname"
split <- parser_args[["split"]] %||% "|"
sep <- parser_args[["sep"]] %||% ","
cores <- as.numeric(parser_args[["cores"]] %||% "1")

# Load libraries
library(LongcellPre)

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load input data
cat("Loading input data...\n")
umi_count <- read.table(umi_count_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_bed <- readRDS(gene_bed_path)

# Load GTF if provided
gtf <- NULL
if (!is.null(gtf_path) && gtf_path != "NO_FILE") {
  cat("Loading GTF annotation...\n")
  gtf <- readRDS(gtf_path)
} else {
  cat("Warning: GTF not provided, skipping isoform imputation.\n")
}

# Run isoform imputation if GTF is available
if (!is.null(gtf)) {
  cat("Running isoform alignment with", cores, "cores...\n")
  UMI_count_to_isoform(
    umi_count = umi_count,
    dir = out_dir,
    gene_bed = gene_bed,
    gtf = gtf,
    gene_col = gene_col,
    bed_gene_col = bed_gene_col,
    bed_strand_col = bed_strand_col,
    filter_only_intron = filter_only_intron,
    mid_offset_thresh = mid_offset_thresh,
    overlap_thresh = overlap_thresh,
    gtf_gene_col = gtf_gene_col,
    gtf_start_col = gtf_start_col,
    gtf_end_col = gtf_end_col,
    gtf_iso_col = gtf_iso_col,
    split = split,
    sep = sep,
    cores = cores
  )
  
  cat("Isoform imputation completed successfully.\n")
} else {
  cat("Skipping isoform imputation (no GTF provided).\n")
}
