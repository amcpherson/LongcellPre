#!/usr/bin/env Rscript

# Split a gene bed RDS into chunks of approximately N genes each.
# Outputs gene_bed_chunk_001.rds, gene_bed_chunk_002.rds, etc.

parse_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    gene_bed_path = NULL,
    genes_per_chunk = 50
  )

  i <- 1
  while(i <= length(args)){
    arg <- args[i]
    if(arg == '--gene_bed_path'){
      opts$gene_bed_path <- args[i+1]
      i <- i + 2
    } else if(arg == '--genes_per_chunk'){
      opts$genes_per_chunk <- as.integer(args[i+1])
      i <- i + 2
    } else {
      stop(paste('Unknown argument:', arg))
    }
  }

  if(is.null(opts$gene_bed_path)) stop('--gene_bed_path is required')

  return(opts)
}

opts <- parse_args()

gene_bed <- readRDS(opts$gene_bed_path)
genes <- unique(gene_bed$gene)
n_genes <- length(genes)
n_chunks <- max(1, ceiling(n_genes / opts$genes_per_chunk))

cat("Splitting", n_genes, "genes into", n_chunks, "chunks for isoform extraction\n")
gene_chunks <- split(genes, ceiling(seq_along(genes) / max(1, n_genes / n_chunks)))

for (i in seq_along(gene_chunks)) {
  chunk_bed <- gene_bed[gene_bed$gene %in% gene_chunks[[i]], ]
  chunk_file <- paste0('gene_bed_chunk_', sprintf('%03d', i), '.rds')
  saveRDS(chunk_bed, chunk_file)
  cat("Wrote chunk", i, "with", length(gene_chunks[[i]]), "genes\n")
}
