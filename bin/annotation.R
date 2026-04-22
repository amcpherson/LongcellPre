#!/usr/bin/env Rscript

library(LongcellPre)
library(GenomicFeatures)
library(tidyr)
library(dplyr)

parse_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(gtf_path = NULL, gene_bed_path = NULL, work_dir = './', overwrite = FALSE)

  i <- 1
  while(i <= length(args)){
    arg <- args[i]
    if(arg == '--gtf_path'){
      opts$gtf_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--gene_bed_path'){
      opts$gene_bed_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--work_dir'){
      opts$work_dir <- args[i+1]
      i <- i + 2
    } else if(arg == '--overwrite'){
      val <- toupper(args[i+1])
      opts$overwrite <- val %in% c('TRUE', 'T', '1')
      i <- i + 2
    } else {
      stop(paste('Unknown argument:', arg))
    }
  }

  if(is.null(opts$gtf_path) && is.null(opts$gene_bed_path)){
    stop('Either --gtf_path or --gene_bed_path must be provided.')
  }

  if(!dir.exists(opts$work_dir)){
    dir.create(opts$work_dir, recursive = TRUE)
  }

  return(opts)
}

opts <- parse_args()

LongcellPre::annotation(gtf_path = opts$gtf_path,
                        gene_bed_path = opts$gene_bed_path,
                        work_dir = opts$work_dir,
                        overwrite = opts$overwrite)
