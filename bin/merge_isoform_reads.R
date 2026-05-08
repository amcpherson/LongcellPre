#!/usr/bin/env Rscript

# Merge per-chunk BarcodeMatchIso files and compute adapter quality metrics.
# Outputs BarcodeMatch/BarcodeMatchIso.txt and BarcodeMatch/adapterNeedle.txt.

library(LongcellPre)
library(dplyr)

parse_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    chunk_dir = '.',
    out_dir = '.'
  )

  i <- 1
  while(i <= length(args)){
    arg <- args[i]
    if(arg == '--chunk_dir'){
      opts$chunk_dir <- args[i+1]
      i <- i + 2
    } else if(arg == '--out_dir'){
      opts$out_dir <- args[i+1]
      i <- i + 2
    } else {
      stop(paste('Unknown argument:', arg))
    }
  }

  return(opts)
}

opts <- parse_args()

dir.create(file.path(opts$out_dir, "BarcodeMatch"), showWarnings = FALSE, recursive = TRUE)

files <- sort(list.files(opts$chunk_dir, pattern = "^chunk_.*\\.txt$", full.names = TRUE))
cat("Merging", length(files), "isoform extraction chunks\n")

chunks <- lapply(files, function(f) {
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})
reads_bc <- as.data.frame(do.call(rbind, chunks))
cat("Combined", nrow(reads_bc), "rows\n")

saveResult(reads_bc, file.path(opts$out_dir, "BarcodeMatch/BarcodeMatchIso.txt"))

if (nrow(reads_bc) > 0) {
  qual <- adapter_dis(data = reads_bc)
  saveResult(qual, file.path(opts$out_dir, "BarcodeMatch/adapterNeedle.txt"))
} else {
  stop("No read is found with valid barcode, please check if your barcode and fastq file match!")
}
cat("Merge completed successfully.\n")
