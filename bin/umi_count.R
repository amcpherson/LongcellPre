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

# Load libraries
library(LongcellPre)

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load input data
cat("Loading input data...\n")
data <- read.table(data_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
qual <- read.table(qual_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
gene_bed <- readRDS(gene_bed_path)

# Get gene strand info
gene_strand <- unique(gene_bed[,c(bed_gene_col,bed_strand_col)])
colnames(gene_strand) <- c("gene", "strand")

# Run UMI count (no internal parallelism - Nextflow handles that)
cat("Running UMI deduplication on chunk...\n")
count <- umi_count(
  cell_exon = data,
  qual = qual,
  gene_strand = gene_strand,
  bar = bar,
  gene = gene,
  isoform = isoform,
  polyA = polyA,
  sim_thresh = sim_thresh,
  split = split,
  sep = sep,
  splice_site_thresh = splice_site_thresh,
  verbose = verbose
)

# Filter and format output
if(!is.null(count) && nrow(count) > 0) {
  count <- as.data.frame(count)
  count <- count[,c("cell","gene","isoform","count","polyA")]
  cat("UMI deduplication completed successfully.\n")
} else {
  # Create empty dataframe with correct column structure
  count <- data.frame(cell = character(),
                      gene = character(),
                      isoform = character(),
                      count = numeric(),
                      polyA = logical(),
                      stringsAsFactors = FALSE)
  cat("Warning: No output generated from UMI count, creating empty table.\n")
}

# Always save output (with headers even if empty)
write.table(count, 
            file = file.path(out_dir, "iso_count.txt"),
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

cat("Output saved to", file.path(out_dir, "iso_count.txt"), "\n")
