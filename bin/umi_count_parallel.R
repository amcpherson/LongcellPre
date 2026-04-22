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
data_path <- parser_args[["data_path"]]
qual_path <- parser_args[["qual_path"]]
gene_bed_path <- parser_args[["gene_bed_path"]]
out_dir <- parser_args[["out_dir"]]

# Extract optional arguments with defaults
bar <- parser_args[["bar"]] %||% "barcode"
gene <- parser_args[["gene"]] %||% "gene"
isoform <- parser_args[["isoform"]] %||% "isoform"
polyA <- parser_args[["polyA"]] %||% "polyA"
sim_thresh <- if (!is.null(parser_args[["sim_thresh"]])) as.numeric(parser_args[["sim_thresh"]]) else NULL
split <- parser_args[["split"]] %||% "|"
sep <- parser_args[["sep"]] %||% ","
splice_site_thresh <- as.numeric(parser_args[["splice_site_thresh"]] %||% "3")
verbose <- parser_args[["verbose"]] %||% "FALSE"
verbose <- as.logical(verbose)
bed_gene_col <- parser_args[["bed_gene_col"]] %||% "gene"
bed_strand_col <- parser_args[["bed_strand_col"]] %||% "strand"
cores <- as.numeric(parser_args[["cores"]] %||% "1")
force_UMI_dedup <- parser_args[["force_UMI_dedup"]] %||% "FALSE"
force_UMI_dedup <- as.logical(force_UMI_dedup)

# Load libraries
library(LongcellPre)

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load input data
cat("Loading input data...\n")
data <- read.table(data_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
qual <- read.table(qual_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_bed <- readRDS(gene_bed_path)

# Run UMI count
cat("Running UMI deduplication with", cores, "cores...\n")
count <- umi_count_parallel(
  data = data,
  qual = qual,
  dir = out_dir,
  gene_bed = gene_bed,
  bar = bar,
  gene = gene,
  isoform = isoform,
  polyA = polyA,
  sim_thresh = sim_thresh,
  split = split,
  sep = sep,
  splice_site_thresh = splice_site_thresh,
  verbose = verbose,
  bed_gene_col = bed_gene_col,
  bed_strand_col = bed_strand_col,
  cores = cores,
  force_UMI_dedup = force_UMI_dedup
)

cat("UMI deduplication completed successfully.\n")
