#!/usr/bin/env Rscript

# Filter gene bed to remove genes without BAM coverage using bedtools.
# Outputs filtered gene_bed as an RDS file.

library(LongcellPre)
library(dplyr)

parse_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    bam_path = NULL,
    gene_bed_path = NULL,
    bedtools = 'bedtools',
    out_path = 'gene_bed_filtered.rds'
  )

  i <- 1
  while(i <= length(args)){
    arg <- args[i]
    if(arg == '--bam_path'){
      opts$bam_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--gene_bed_path'){
      opts$gene_bed_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--bedtools'){
      opts$bedtools <- args[i+1]
      i <- i + 2
    } else if(arg == '--out_path'){
      opts$out_path <- args[i+1]
      i <- i + 2
    } else {
      stop(paste('Unknown argument:', arg))
    }
  }

  if(is.null(opts$bam_path)) stop('--bam_path is required')
  if(is.null(opts$gene_bed_path)) stop('--gene_bed_path is required')

  return(opts)
}

opts <- parse_args()

gene_bed <- readRDS(opts$gene_bed_path)

gene_range <- gene_bed %>%
  group_by(gene) %>%
  summarise(chr = unique(chr), start = min(start), end = max(end), strand = unique(strand))
gene_range <- gene_range[, c("chr", "start", "end", "strand", "gene")]
write.table(gene_range, "gene_range.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

noncover <- bamGeneCoverage(
  bam = opts$bam_path,
  gene_range_bed = "gene_range.txt",
  outdir = ".",
  bedtools = opts$bedtools
)
if (!is.null(noncover)) {
  gene_bed <- gene_bed %>% filter(!gene %in% noncover$gene)
}
cat("Filtered gene bed:", length(unique(gene_bed$gene)), "genes remaining\n")
saveRDS(gene_bed, opts$out_path)
