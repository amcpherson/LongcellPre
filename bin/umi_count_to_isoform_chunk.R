#!/usr/bin/env Rscript

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

parser_args <- list()
for (i in seq(1, length(args), by = 2)) {
  if (i < length(args)) {
    param_name <- gsub("^--", "", args[i])
    param_value <- args[i + 1]
    parser_args[[param_name]] <- param_value
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Extract arguments
umi_count_path     <- parser_args[["umi_count_path"]]
gtf_path           <- parser_args[["gtf_path"]]
out_path           <- parser_args[["out_path"]]
filter_only_intron <- as.logical(parser_args[["filter_only_intron"]] %||% "TRUE")
mid_offset_thresh  <- as.numeric(parser_args[["mid_offset_thresh"]]  %||% "3")
overlap_thresh     <- as.numeric(parser_args[["overlap_thresh"]]     %||% "0")
gtf_gene_col       <- parser_args[["gtf_gene_col"]]  %||% "gene"
gtf_start_col      <- parser_args[["gtf_start_col"]] %||% "start"
gtf_end_col        <- parser_args[["gtf_end_col"]]   %||% "end"
gtf_iso_col        <- parser_args[["gtf_iso_col"]]   %||% "transname"
split              <- parser_args[["split"]] %||% "|"
sep                <- parser_args[["sep"]]   %||% ","

library(LongcellPre)

cat("Loading input data...\n")
umi_count <- read.table(umi_count_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Chunk has", nrow(umi_count), "rows across", length(unique(umi_count$gene)), "genes\n")

cat("Loading GTF annotation...\n")
gtf <- readRDS(gtf_path)

cat("Running isoform alignment...\n")
result <- cells_genes_isos_count(
  data               = umi_count,
  gtf                = gtf,
  thresh             = mid_offset_thresh,
  overlap_thresh     = overlap_thresh,
  filter_only_intron = filter_only_intron,
  gtf_gene_col       = gtf_gene_col,
  gtf_iso_col        = gtf_iso_col,
  gtf_start_col      = gtf_start_col,
  gtf_end_col        = gtf_end_col,
  split              = split,
  sep                = sep
)

if (!is.null(result) && nrow(result) > 0) {
  write.table(result, file = out_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  cat("Written", nrow(result), "rows to", out_path, "\n")
} else {
  # Write empty file with header so merge step doesn't fail
  write.table(
    data.frame(cell = character(), isoform = character(), count = numeric(), gene = character()),
    file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  cat("No results for this chunk, wrote empty file\n")
}
