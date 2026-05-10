#!/usr/bin/env Rscript

# Generate UMI-collapsed FASTQ for a chunk of genes (no mapping).
# Called per-chunk by Nextflow; the merged FASTQ is mapped in a separate step.

args <- commandArgs(trailingOnly = TRUE)

parser_args <- list()
for (i in seq(1, length(args), by = 2)) {
  if (i < length(args)) {
    param_name <- gsub("^--", "", args[i])
    param_value <- args[i + 1]
    parser_args[[param_name]] <- param_value
  }
}

umi_count_path <- parser_args[["umi_count_path"]]
gene_bed_path  <- parser_args[["gene_bed_path"]]
genome_name    <- parser_args[["genome_name"]]
out_fastq      <- parser_args[["out_fastq"]]
out_annot      <- parser_args[["out_annot"]]

library(LongcellPre)

options(future.globals.maxSize = 24 * 1024^3)

umi_count <- read.table(umi_count_path, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
gene_bed  <- readRDS(gene_bed_path)

cat("Loading genome...\n")
genome <- load_genome(genome_name)

cat("Generating FASTQ for", length(unique(umi_count$gene)), "genes...\n")
reads <- isoformCount2Reads(umi_count, genome, gene_bed, out_fastq)

qname <- as.character(reads@id)
annot <- extractAnnotFromQname(qname, "cell")
qname <- cbind(qname, annot)

write.table(qname, file = out_annot, sep = ",", quote = FALSE,
            row.names = FALSE)

cat("Chunk FASTQ generation complete.\n")
