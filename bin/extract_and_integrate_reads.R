#!/usr/bin/env Rscript

# This script performs the isoform extraction and integration step only.
# Barcode extraction and mapping are handled by earlier Nextflow processes.
# It mirrors the isoform extraction section of reads_extract_bc() in pipeline.R.

library(LongcellPre)
library(dplyr)
library(future)
library(future.apply)

parse_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    barcode_path = NULL,
    bam_path = NULL,
    gene_bed_path = NULL,
    genome_name = NULL,
    work_dir = './',
    toolkit = 5,
    bedtools = 'bedtools',
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
    } else if(arg == '--bedtools'){
      opts$bedtools <- args[i+1]
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

# Set up parallel plan for reads_extraction (uses future_lapply internally)
cores <- LongcellPre::coreDetect(opts$cores)
if(cores > 1){
  plan(strategy = "multisession", workers = cores)
} else {
  plan(strategy = "sequential")
}

# Create output directories
dir.create(file.path(opts$work_dir, "BarcodeMatch"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opts$work_dir, "annotation"), showWarnings = FALSE, recursive = TRUE)

# Load inputs
cat("Loading gene bed annotation...\n")
gene_bed <- readRDS(opts$gene_bed_path)

cat("Loading barcode match data...\n")
bc <- read.table(opts$barcode_path, header = TRUE, sep = "\t")

# Filter genes without coverage (bamGeneCoverage)
cat("Filtering genes without read coverage...\n")
gene_range <- gene_bed %>%
  group_by(gene) %>%
  summarise(chr = unique(chr), start = min(start), end = max(end), strand = unique(strand))
gene_range <- gene_range[, c("chr", "start", "end", "strand", "gene")]
gene_range_file <- file.path(opts$work_dir, "annotation/gene_range.txt")
write.table(gene_range, gene_range_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

noncover <- bamGeneCoverage(
  bam = opts$bam_path,
  gene_range_bed = gene_range_file,
  outdir = file.path(opts$work_dir, "annotation"),
  bedtools = opts$bedtools
)
if(!is.null(noncover)){
  gene_bed <- gene_bed %>% filter(!gene %in% noncover$gene)
}

# Extract isoform information from BAM
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

saveResult(reads_bc, file.path(opts$work_dir, "BarcodeMatch/BarcodeMatchIso.txt"))

if(nrow(reads_bc) > 0){
  # Evaluate data quality
  qual <- adapter_dis(data = reads_bc)
  saveResult(qual, file.path(opts$work_dir, "BarcodeMatch/adapterNeedle.txt"))
} else {
  stop("No read is found with valid barcode, please check if your barcode and fastq file match!")
}

cat("Isoform extraction and integration completed successfully.\n")
