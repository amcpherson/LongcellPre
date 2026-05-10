#!/usr/bin/env Rscript

# This script performs the isoform extraction and integration step for a
# gene chunk. Gene coverage filtering and gene-bed splitting are handled
# by earlier Nextflow processes (FILTER_GENE_BED / SPLIT_GENE_BED).
# Each chunk is processed sequentially — Nextflow handles parallelism
# across chunks.

library(LongcellPre)
library(dplyr)
library(future)
library(future.apply)

# Increase future globals size limit for parallel workers
options(future.globals.maxSize = 24 * 1024^3)

parse_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    barcode_path = NULL,
    bam_path = NULL,
    gene_bed_path = NULL,
    genome_name = NULL,
    work_dir = './',
    toolkit = 5,
    map_qual = 30,
    end_flank = 200,
    splice_site_bin = 2,
    cores = 1
  )

  i <- 1
  while(i <= length(args)){
    arg <- args[i]
    if(arg == '--barcode_path'){
      opts$barcode_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--bam_path'){
      opts$bam_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--gene_bed_path'){
      opts$gene_bed_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--genome_name'){
      opts$genome_name <- args[i+1]
      i <- i + 2
    } else if(arg == '--work_dir'){
      opts$work_dir <- args[i+1]
      i <- i + 2
    } else if(arg == '--toolkit'){
      opts$toolkit <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--map_qual'){
      opts$map_qual <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--end_flank'){
      opts$end_flank <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--splice_site_bin'){
      opts$splice_site_bin <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--cores'){
      opts$cores <- as.integer(args[i+1])
      i <- i + 2
    } else {
      stop(paste('Unknown argument:', arg))
    }
  }

  if(is.null(opts$barcode_path)){
    stop('--barcode_path is required')
  }
  if(is.null(opts$bam_path)){
    stop('--bam_path is required')
  }
  if(is.null(opts$gene_bed_path)){
    stop('--gene_bed_path is required')
  }
  if(is.null(opts$genome_name)){
    stop('--genome_name is required')
  }

  return(opts)
}

opts <- parse_args()

# Use sequential plan — Nextflow handles parallelism across gene chunks.
# Only use multisession if cores > 1 is explicitly requested for within-chunk
# parallelism (e.g. when running outside Nextflow).
cores <- LongcellPre::coreDetect(opts$cores)
if(cores > 1){
  plan(strategy = "multisession", workers = cores)
} else {
  plan(strategy = "sequential")
}

# Load inputs
cat("Loading gene bed chunk...\n")
gene_bed <- readRDS(opts$gene_bed_path)
cat("Gene chunk contains", length(unique(gene_bed$gene)), "genes\n")

cat("Loading barcode match data...\n")
bc <- read.table(opts$barcode_path, header = TRUE, sep = "\t")

# Extract isoform information from BAM for this gene chunk
cat("Extracting isoforms from BAM...\n")
suppressWarnings({genome <- load_genome(opts$genome_name)})
reads <- reads_extraction(
  bam_path = opts$bam_path,
  gene_bed = gene_bed,
  genome = genome,
  toolkit = opts$toolkit,
  map_qual = opts$map_qual,
  end_flank = opts$end_flank,
  splice_site_bin = opts$splice_site_bin
)

# Integrate barcode data with isoform data
cat("Integrating barcode and isoform data...\n")
reads_bc <- inner_join(bc, reads, by = c("name" = "qname"))
reads_bc <- reads_bc %>%
  mutate(polyA.x = as.numeric(polyA.x), polyA.y = as.numeric(polyA.y)) %>%
  mutate(polyA = polyA.x & polyA.y) %>%
  dplyr::select(-polyA.x, -polyA.y)

out_path <- file.path(opts$work_dir, "BarcodeMatchIso_chunk.txt")
saveResult(reads_bc, out_path)
cat("Wrote", nrow(reads_bc), "rows to", out_path, "\n")
cat("Chunk extraction and integration completed successfully.\n")
