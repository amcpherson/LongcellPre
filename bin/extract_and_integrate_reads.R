#!/usr/bin/env Rscript

library(LongcellPre)

parse_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    fastq_path = NULL,
    barcode_path = NULL,
    annotation_dir = NULL,
    gene_bed_path = NULL,
    work_dir = './',
    genome_path = NULL,
    genome_name = NULL,
    toolkit = 5,
    protocol = '10X',
    adapter = NULL,
    minimap_bed_path = NULL,
    window = 10,
    step = 2,
    left_flank = 0,
    right_flank = 0,
    drop_adapter = FALSE,
    polyA_bin = 20,
    polyA_base_count = 15,
    polyA_len = 10,
    barcode_len = 16,
    mu = 15,
    sigma = 10,
    k = 6,
    batch = 100,
    top = 5,
    cos_thresh = 0.25,
    alpha = 0.05,
    edit_thresh = 3,
    mean_edit_thresh = 1.5,
    UMI_len = 10,
    UMI_flank = 1,
    minimap2 = 'minimap2',
    samtools = 'samtools',
    bedtools = 'bedtools',
    map_qual = 30,
    end_flank = 200,
    splice_site_bin = 2,
    force_barcode_match = FALSE,
    force_map = FALSE,
    force_isoform_extract = FALSE,
    force_rerun = FALSE,
    cores = 1
  )

  i <- 1
  while(i <= length(args)){
    arg <- args[i]
    if(arg == '--fastq_path'){
      opts$fastq_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--barcode_path'){
      opts$barcode_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--annotation_dir'){
      opts$annotation_dir <- args[i+1]
      i <- i + 2
    } else if(arg == '--gene_bed_path'){
      opts$gene_bed_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--work_dir'){
      opts$work_dir <- args[i+1]
      i <- i + 2
    } else if(arg == '--genome_path'){
      opts$genome_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--genome_name'){
      opts$genome_name <- args[i+1]
      i <- i + 2
    } else if(arg == '--toolkit'){
      opts$toolkit <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--protocol'){
      opts$protocol <- args[i+1]
      i <- i + 2
    } else if(arg == '--adapter'){
      opts$adapter <- args[i+1]
      i <- i + 2
    } else if(arg == '--minimap_bed_path'){
      opts$minimap_bed_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--window'){
      opts$window <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--step'){
      opts$step <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--left_flank'){
      opts$left_flank <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--right_flank'){
      opts$right_flank <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--drop_adapter'){
      val <- toupper(args[i+1])
      opts$drop_adapter <- val %in% c('TRUE', 'T', '1')
      i <- i + 2
    } else if(arg == '--polyA_bin'){
      opts$polyA_bin <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--polyA_base_count'){
      opts$polyA_base_count <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--polyA_len'){
      opts$polyA_len <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--barcode_len'){
      opts$barcode_len <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--mu'){
      opts$mu <- as.numeric(args[i+1])
      i <- i + 2
    } else if(arg == '--sigma'){
      opts$sigma <- as.numeric(args[i+1])
      i <- i + 2
    } else if(arg == '--k'){
      opts$k <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--batch'){
      opts$batch <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--top'){
      opts$top <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--cos_thresh'){
      opts$cos_thresh <- as.numeric(args[i+1])
      i <- i + 2
    } else if(arg == '--alpha'){
      opts$alpha <- as.numeric(args[i+1])
      i <- i + 2
    } else if(arg == '--edit_thresh'){
      opts$edit_thresh <- as.numeric(args[i+1])
      i <- i + 2
    } else if(arg == '--mean_edit_thresh'){
      opts$mean_edit_thresh <- as.numeric(args[i+1])
      i <- i + 2
    } else if(arg == '--UMI_len'){
      opts$UMI_len <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--UMI_flank'){
      opts$UMI_flank <- as.integer(args[i+1])
      i <- i + 2
    } else if(arg == '--minimap2'){
      opts$minimap2 <- args[i+1]
      i <- i + 2
    } else if(arg == '--samtools'){
      opts$samtools <- args[i+1]
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
    } else if(arg == '--force_barcode_match'){
      val <- toupper(args[i+1])
      opts$force_barcode_match <- val %in% c('TRUE', 'T', '1')
      i <- i + 2
    } else if(arg == '--force_map'){
      val <- toupper(args[i+1])
      opts$force_map <- val %in% c('TRUE', 'T', '1')
      i <- i + 2
    } else if(arg == '--force_isoform_extract'){
      val <- toupper(args[i+1])
      opts$force_isoform_extract <- val %in% c('TRUE', 'T', '1')
      i <- i + 2
    } else if(arg == '--force_rerun'){
      val <- toupper(args[i+1])
      opts$force_rerun <- val %in% c('TRUE', 'T', '1')
      i <- i + 2
    } else if(arg == '--cores'){
      opts$cores <- as.integer(args[i+1])
      i <- i + 2
    } else {
      stop(paste('Unknown argument:', arg))
    }
  }

  if(is.null(opts$fastq_path)){
    stop('--fastq_path is required')
  }
  if(is.null(opts$barcode_path)){
    stop('--barcode_path is required')
  }
  if(is.null(opts$genome_path)){
    stop('--genome_path is required')
  }
  if(is.null(opts$genome_name)){
    stop('--genome_name is required')
  }
  if(is.null(opts$annotation_dir) && is.null(opts$gene_bed_path)){
    stop('Either --annotation_dir or --gene_bed_path must be provided')
  }

  if(!dir.exists(opts$work_dir)){
    dir.create(opts$work_dir, recursive = TRUE)
  }

  return(opts)
}

opts <- parse_args()

if(!is.null(opts$annotation_dir)){
  gene_bed <- readRDS(file.path(opts$annotation_dir, 'gene_bed.rds'))
} else {
  gene_bed <- read.table(opts$gene_bed_path, header = TRUE)
  if(!all(c('chr','start','end','strand','gene') %in% colnames(gene_bed))){
    stop('gene_bed_path must contain chr, start, end, strand, gene columns')
  }
}

LongcellPre::reads_extract_bc(
  fastq_path = opts$fastq_path,
  barcode_path = opts$barcode_path,
  gene_bed = gene_bed,
  adapter = opts$adapter,
  genome_path = opts$genome_path,
  genome_name = opts$genome_name,
  toolkit = opts$toolkit,
  protocol = opts$protocol,
  minimap_bed_path = opts$minimap_bed_path,
  work_dir = opts$work_dir,
  window = opts$window,
  step = opts$step,
  left_flank = opts$left_flank,
  right_flank = opts$right_flank,
  drop_adapter = opts$drop_adapter,
  polyA_bin = opts$polyA_bin,
  polyA_base_count = opts$polyA_base_count,
  polyA_len = opts$polyA_len,
  barcode_len = opts$barcode_len,
  mu = opts$mu,
  sigma = opts$sigma,
  k = opts$k,
  batch = opts$batch,
  top = opts$top,
  cos_thresh = opts$cos_thresh,
  alpha = opts$alpha,
  edit_thresh = opts$edit_thresh,
  mean_edit_thresh = opts$mean_edit_thresh,
  UMI_len = opts$UMI_len,
  UMI_flank = opts$UMI_flank,
  minimap2 = opts$minimap2,
  samtools = opts$samtools,
  bedtools = opts$bedtools,
  map_qual = opts$map_qual,
  end_flank = opts$end_flank,
  splice_site_bin = opts$splice_site_bin,
  force_barcode_match = opts$force_barcode_match,
  force_map = opts$force_map,
  force_isoform_extract = opts$force_isoform_extract,
  force_rerun = opts$force_rerun,
  cores = opts$cores
)
