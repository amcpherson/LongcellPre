#!/usr/bin/env Rscript

library(LongcellPre)

parse_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    fastq_path = NULL,
    barcode_path = NULL,
    out_name = "polish.fq.gz",
    toolkit = 5,
    protocol = "10X",
    adapter = NULL,
    window = 10,
    step = 2,
    left_flank = 0,
    right_flank = 0,
    drop_adapter = FALSE,
    polyA_bin = 20,
    polyA_base_count = 15,
    polyA_len = 10,
    barcode_len = 16,
    mu = 20,
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
    } else if(arg == '--out_name'){
      opts$out_name <- args[i+1]
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

  return(opts)
}

opts <- parse_args()

# Call the extractTagBc function with parsed arguments
bc <- LongcellPre::extractTagBc(
  fastq_path = opts$fastq_path,
  barcode_path = opts$barcode_path,
  out_name = opts$out_name,
  toolkit = opts$toolkit,
  protocol = opts$protocol,
  adapter = opts$adapter,
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
  cores = opts$cores
)

# Save the barcode data
write.table(bc, "BarcodeMatch.txt", sep = "\t", quote = FALSE, row.names = FALSE)
