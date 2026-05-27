nextflow.enable.dsl=2

// ============================================================================
// LongcellPre Nextflow Pipeline
// ============================================================================
// This pipeline processes long-read single-cell RNA-seq data with barcode
// and UMI tagging.
//
// FASTQ Chunking:
// The input FASTQ is split into chunks of ~1M reads each (configurable via
// params.fastq_chunk_size). Each chunk is processed independently for
// barcode/UMI extraction, then results are merged. This reduces peak memory
// usage without affecting accuracy.
//
// UMI Counting Chunking:
// The barcode/isoform data is split into chunks balanced by multi-exon read
// count (~50k multi-exon reads per chunk, configurable via
// params.multiexon_reads_per_chunk). Genes are assigned via greedy bin-packing
// so that each chunk has approximately equal workload. Each gene chunk is
// processed independently for UMI deduplication, then results are merged.
//
// ============================================================================

params.gtf_path = null
params.gene_bed_path = null
params.fastq_path = null
params.barcode_path = null
params.genome_path = null
params.genome_name = null
params.minimap_bed_path = null
params.toolkit = 5
params.protocol = '10X'
params.bedtools = 'bedtools'
params.map_qual = 30
params.end_flank = 200
params.splice_site_bin = 2
params.mean_edit_thresh = 1.5
params.results_dir = 'results'
params.overwrite = true
params.cores = 4
params.to_isoform = true
// Barcode extraction parameters
params.adapter = null
params.window = 10
params.step = 2
params.left_flank = 0
params.right_flank = 0
params.drop_adapter = false
params.polyA_bin = 20
params.polyA_base_count = 15
params.polyA_len = 10
params.barcode_len = 16
params.mu = 15
params.sigma = 10
params.k = 6
params.batch = 100
params.top = 5
params.cos_thresh = 0.25
params.alpha = 0.05
params.edit_thresh = 3
params.UMI_len = 10
params.UMI_flank = 1
// UMI deduplication parameters
params.splice_site_thresh = 3
params.sim_thresh = null
params.verbose = false
params.bed_gene_col = "gene"
params.bed_strand_col = "strand"
// Isoform imputation parameters
params.filter_only_intron = true
params.mid_offset_thresh = 3
params.overlap_thresh = 0
params.gtf_gene_col = "gene"
params.gtf_start_col = "start"
params.gtf_end_col = "end"
params.gtf_iso_col = "transname"
params.split = "|"
params.sep = ","
// FASTQ chunking parameters
params.fastq_chunk_size = 1000000  // Number of reads per chunk (1M default)
// UMI counting chunking parameters
params.multiexon_reads_per_chunk = 50000  // Target multi-exon reads per chunk
// Tool paths
params.minimap2 = "minimap2"
params.samtools = "samtools"

process RUN_ANNOTATION {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path gtf_file
    path gene_bed_file

    output:
    path 'annotation/gene_bed.rds'
    path 'annotation/exon_gtf.rds'

    script:
    def gtf_arg = gtf_file.name != 'NO_FILE' ? "--gtf_path ${gtf_file}" : ''
    def gene_bed_arg = gene_bed_file.name != 'NO_FILE' ? "--gene_bed_path ${gene_bed_file}" : ''
    def overwrite_arg = params.overwrite ? 'TRUE' : 'FALSE'
    """
    mkdir -p annotation
    annotation.R ${gtf_arg} ${gene_bed_arg} --work_dir . --overwrite ${overwrite_arg}
    """
}

process SPLIT_FASTQ {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path fastq_file

    output:
    path 'chunk_*.fastq.gz'

    script:
    """
    # Split FASTQ into chunks (4 lines per FASTQ record)
    gunzip -c '${fastq_file}' | split -l \$((${params.fastq_chunk_size} * 4)) -a 4 --numeric-suffixes=1 -d - chunk_raw_
    
    # Compress each chunk
    for f in chunk_raw_*; do
        if [ -f "\$f" ]; then
            cat "\$f" | gzip > chunk_\${f##chunk_raw_}.fastq.gz
            rm "\$f"
        fi
    done
    
    # Verify chunks were created
    if ! ls chunk_*.fastq.gz 1> /dev/null 2>&1; then
        echo "Error: No chunks created. Check input FASTQ file."
        exit 1
    fi
    """
}

process EXTRACT_TAG_BC_CHUNK {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path fastq_chunk
    path barcode_file

    output:
    path 'polish.fq.gz'
    path 'BarcodeMatch.txt'

    script:
    def adapter_arg = params.adapter ? "--adapter ${params.adapter}" : ''
    """
    extract_tag_bc.R \
      --fastq_path ${fastq_chunk} \
      --barcode_path ${barcode_file} \
      --toolkit ${params.toolkit} \
      --protocol ${params.protocol} \
      ${adapter_arg} \
      --window ${params.window} \
      --step ${params.step} \
      --left_flank ${params.left_flank} \
      --right_flank ${params.right_flank} \
      --drop_adapter ${params.drop_adapter} \
      --polyA_bin ${params.polyA_bin} \
      --polyA_base_count ${params.polyA_base_count} \
      --polyA_len ${params.polyA_len} \
      --barcode_len ${params.barcode_len} \
      --mu ${params.mu} \
      --sigma ${params.sigma} \
      --k ${params.k} \
      --batch ${params.batch} \
      --top ${params.top} \
      --cos_thresh ${params.cos_thresh} \
      --alpha ${params.alpha} \
      --edit_thresh ${params.edit_thresh} \
      --mean_edit_thresh ${params.mean_edit_thresh} \
      --UMI_len ${params.UMI_len} \
      --UMI_flank ${params.UMI_flank} \
      --cores ${params.cores}
    """
}

process MERGE_BARCODE_MATCH {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path 'barcode_*.txt'

    output:
    path 'BarcodeMatch.txt'

    publishDir "${params.results_dir}", mode: 'copy', overwrite: true

    script:
    """
    # Merge barcode match files, keeping header from first file
    awk 'FNR==1 && NR>1 { next; } { print }' barcode_*.txt > BarcodeMatch.txt
    """
}

process MAP_POLISHED_FASTQ_CHUNK {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path fastq_file
    path genome_file
    path bed_file

    output:
    path '*.bam'

    script:
    def bed_arg = bed_file.name != 'NO_FILE' ? "--junc-bed ${bed_file}" : ''
    """
    ${params.minimap2} -ax splice -uf --sam-hit-only -t ${params.cores} ${genome_file} ${fastq_file} ${bed_arg} | \
    ${params.samtools} view -bS -@ ${params.cores} - | \
    ${params.samtools} sort - -@ ${params.cores} -o chunk.bam
    """
}

process MERGE_BAM {
    publishDir "${params.results_dir}/bam", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path 'chunk_*.bam'

    output:
    path 'polish.bam'
    path 'polish.bam.bai'

    script:
    """
    ${params.samtools} merge -@ ${params.cores} polish.bam chunk_*.bam
    ${params.samtools} index polish.bam
    """
}

process FILTER_GENE_BED {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path bam_file
    path bai_file
    path gene_bed_rds

    output:
    path 'gene_bed_filtered.rds'

    script:
    """
    filter_gene_bed.R \
      --bam_path ${bam_file} \
      --gene_bed_path ${gene_bed_rds} \
      --bedtools ${params.bedtools} \
      --out_path gene_bed_filtered.rds
    """
}


process EXTRACT_AND_INTEGRATE_READS_CHUNK {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path barcode_file
    path bam_file
    path gene_bed_chunk

    output:
    path 'BarcodeMatchIso_chunk.txt'

    script:
    """
    ${params.samtools} index ${bam_file}

    extract_and_integrate_reads.R \
      --barcode_path ${barcode_file} \
      --bam_path ${bam_file} \
      --gene_bed_path ${gene_bed_chunk} \
      --genome_name ${params.genome_name} \
      --work_dir . \
      --toolkit ${params.toolkit} \
      --map_qual ${params.map_qual} \
      --end_flank ${params.end_flank} \
      --splice_site_bin ${params.splice_site_bin} \
      --cores 1
    """
}

process MERGE_ISOFORM_READS {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path 'chunk_*.txt'

    output:
    path 'BarcodeMatch/BarcodeMatchIso.txt'
    path 'BarcodeMatch/adapterNeedle.txt'

    script:
    """
    merge_isoform_reads.R \
      --chunk_dir . \
      --out_dir .
    """
}

process SPLIT_UMI_DATA {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path barcode_iso_file
    path adapter_needle_file

    output:
    path 'data_chunk_*.txt'
    path 'qual_chunk_*.txt'

    script:
    """
    #!/usr/bin/env Rscript
    
    # Load data files
    data <- read.table('${barcode_iso_file}', header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
    qual <- read.table('${adapter_needle_file}', header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
    
    # Count multi-exon reads per gene (these drive runtime/memory)
    is_multi <- grepl("|", data\$isoform, fixed = TRUE)
    gene_multi_counts <- tapply(is_multi, data\$gene, sum)
    gene_multi_counts <- sort(gene_multi_counts, decreasing = TRUE)
    genes <- names(gene_multi_counts)
    
    n_genes <- length(genes)
    n_chunks <- max(1, ceiling(sum(gene_multi_counts) / ${params.multiexon_reads_per_chunk}))
    
    cat("Splitting", n_genes, "genes into", n_chunks, "chunks (balanced by multi-exon reads)\\n")
    cat("Total multi-exon reads:", sum(gene_multi_counts), "\\n")
    cat("Target per chunk:", ${params.multiexon_reads_per_chunk}, "\\n")
    cat("Top 5 genes by multi-exon reads:\\n")
    for (g in head(genes, 5)) {
        cat("  ", g, ":", gene_multi_counts[g], "\\n")
    }
    
    # Greedy bin-packing: assign each gene (largest first) to the least-loaded chunk
    chunk_totals <- rep(0, n_chunks)
    gene_assignments <- integer(n_genes)
    for (i in seq_along(genes)) {
        min_chunk <- which.min(chunk_totals)
        gene_assignments[i] <- min_chunk
        chunk_totals[min_chunk] <- chunk_totals[min_chunk] + gene_multi_counts[i]
    }
    
    # Write each chunk
    for (i in seq_len(n_chunks)) {
        chunk_genes <- genes[gene_assignments == i]
        chunk_data <- data[data\$gene %in% chunk_genes, ]
        
        chunk_file <- paste0('data_chunk_', sprintf('%03d', i), '.txt')
        write.table(chunk_data, file = chunk_file, sep = "\t", quote = FALSE, 
                    row.names = FALSE, col.names = TRUE)
        cat("Chunk", i, ":", length(chunk_genes), "genes,", nrow(chunk_data), "reads,",
            chunk_totals[i], "multi-exon\\n")
    }
    
    # Qual file is the same for all chunks (no need to split), but copy for each chunk
    for (i in seq_len(n_chunks)) {
        qual_file <- paste0('qual_chunk_', sprintf('%03d', i), '.txt')
        write.table(qual, file = qual_file, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
    }
    
    cat("Done splitting UMI data into", n_chunks, "chunks\\n")
    """
}

process UMI_COUNT_PARALLEL_CHUNK {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path barcode_iso_chunk
    path adapter_needle_chunk
    path gene_bed_file

    output:
    path 'iso_count_*.txt'

    script:
    def sim_thresh_arg = params.sim_thresh ? "--sim_thresh ${params.sim_thresh}" : ''
    def verbose_arg = params.verbose ? "--verbose TRUE" : "--verbose FALSE"
    """
    # Extract chunk number from filename
    chunk_num=\$(basename '${barcode_iso_chunk}' | sed 's/data_chunk_\\([0-9]*\\).txt/\\1/')
    
    # Process this gene chunk with simplified umi_count.R (no internal parallelism)
    # Nextflow handles parallelism at chunk level
    umi_count.R \
      --data_path ${barcode_iso_chunk} \
      --qual_path ${adapter_needle_chunk} \
      --gene_bed_path ${gene_bed_file} \
      --out_dir . \
      --splice_site_thresh ${params.splice_site_thresh} \
      ${sim_thresh_arg} \
      ${verbose_arg} \
      --bed_gene_col ${params.bed_gene_col} \
      --bed_strand_col ${params.bed_strand_col}
    
    # Rename output to include chunk number
    if [ -f iso_count.txt ]; then
        mv iso_count.txt iso_count_\${chunk_num}.txt
    fi
    """
}

process MERGE_UMI_COUNT {
    publishDir "${params.results_dir}/out", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path 'iso_count_*.txt'

    output:
    path 'iso_count.txt'

    script:
    """
    # Merge all iso_count chunks, keeping header from first file
    awk 'FNR==1 && NR>1 { next; } { print }' iso_count_*.txt > iso_count.txt
    
    # Verify output
    lines=\$(wc -l < iso_count.txt)
    echo "Merged UMI count file with \$lines lines"
    """
}

process SPLIT_ISO_INPUT {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path umi_count_file

    output:
    path 'iso_chunk_*.txt'

    script:
    """
    #!/usr/bin/env Rscript

    data <- read.table('${umi_count_file}', header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

    genes <- unique(data\$gene)
    n_genes <- length(genes)
    n_chunks <- max(1, ceiling(n_genes / 500))

    cat("Splitting", n_genes, "genes into", n_chunks, "chunks\\n")

    gene_chunks <- split(genes, ceiling(seq_along(genes) / max(1, n_genes / n_chunks)))

    for (i in seq_along(gene_chunks)) {
        chunk_data <- data[data\$gene %in% gene_chunks[[i]], ]
        chunk_file <- paste0('iso_chunk_', sprintf('%03d', i), '.txt')
        write.table(chunk_data, file = chunk_file, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
        cat("Wrote", nrow(chunk_data), "rows to", chunk_file, "\\n")
    }
    """
}

process UMI_COUNT_TO_ISOFORM_CHUNK {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path umi_count_chunk
    path gtf_file

    output:
    path 'iso_count_mat_chunk.txt'

    script:
    def filter_arg = params.filter_only_intron ? "--filter_only_intron TRUE" : "--filter_only_intron FALSE"
    """
    umi_count_to_isoform_chunk.R \
      --umi_count_path ${umi_count_chunk} \
      --gtf_path ${gtf_file} \
      --out_path iso_count_mat_chunk.txt \
      ${filter_arg} \
      --mid_offset_thresh ${params.mid_offset_thresh} \
      --overlap_thresh ${params.overlap_thresh} \
      --gtf_gene_col ${params.gtf_gene_col} \
      --gtf_start_col ${params.gtf_start_col} \
      --gtf_end_col ${params.gtf_end_col} \
      --gtf_iso_col ${params.gtf_iso_col} \
      --split '${params.split}' \
      --sep '${params.sep}'
    """
}

process MERGE_ISO_MAT {
    publishDir "${params.results_dir}/out", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path 'iso_count_mat_*.txt'

    output:
    path 'gene'
    path 'isoform'

    script:
    """
    #!/usr/bin/env Rscript

    library(LongcellPre)

    files <- list.files(pattern = "^iso_count_mat_.*\\\\.txt\$")
    cat("Merging", length(files), "isoform count chunks\\n")

    chunks <- lapply(files, function(f) {
        read.table(f, header = TRUE, sep = "\\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
    })
    combined <- as.data.frame(do.call(rbind, chunks))
    # Drop empty-header rows from chunks with no results
    combined <- combined[nchar(combined\$cell) > 0, ]
    cat("Combined", nrow(combined), "rows\\n")

    saveIsoMat(combined, ".")
    cat("Isoform matrix saved\\n")
    """
}

process SPLIT_CONSENSUS_INPUT {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path umi_count_file

    output:
    path 'consensus_chunk_*.txt'

    script:
    """
    #!/usr/bin/env Rscript

    data <- read.table('${umi_count_file}', header = TRUE, sep = "\\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

    genes <- unique(data\$gene)
    n_genes <- length(genes)
    n_chunks <- max(1, ceiling(n_genes / 500))

    cat("Splitting", n_genes, "genes into", n_chunks, "consensus chunks\\n")

    gene_chunks <- split(genes, ceiling(seq_along(genes) / max(1, n_genes / n_chunks)))

    for (i in seq_along(gene_chunks)) {
        chunk_data <- data[data\$gene %in% gene_chunks[[i]], ]
        chunk_file <- paste0('consensus_chunk_', sprintf('%03d', i), '.txt')
        write.table(chunk_data, file = chunk_file, sep = "\\t", quote = FALSE,
                    row.names = FALSE, col.names = TRUE)
        cat("Wrote", nrow(chunk_data), "rows to", chunk_file, "\\n")
    }
    """
}

process UMI_CONSENSUS_FASTQ_CHUNK {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path umi_count_chunk
    path gene_bed_file

    output:
    path 'chunk.fq.gz'
    path 'chunk_annot.csv'

    script:
    """
    umi_consensus_fastq_chunk.R \
      --umi_count_path ${umi_count_chunk} \
      --gene_bed_path ${gene_bed_file} \
      --genome_name ${params.genome_name} \
      --out_fastq chunk.fq.gz \
      --out_annot chunk_annot.csv
    """
}

process MERGE_CONSENSUS_AND_MAP {
    publishDir "${params.results_dir}/out", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path 'chunk_*.fq.gz'
    path 'annot_*.csv'
    path genome_file
    path bed_file

    output:
    path 'UMI_collapsed.fq.gz'
    path 'reads_annot.csv'
    path 'UMI_collapsed.bam'
    path 'UMI_collapsed.bam.bai'

    script:
    def bed_arg = bed_file.name != 'NO_FILE' ? "--junc-bed ${bed_file}" : ''
    """
    # Merge FASTQ chunks
    cat chunk_*.fq.gz > UMI_collapsed.fq.gz

    # Merge annotation CSVs (keep header from first file only)
    awk 'FNR==1 && NR>1 { next; } { print }' annot_*.csv > reads_annot.csv

    # Map merged FASTQ to genome
    ${params.minimap2} -ax splice -uf --sam-hit-only -t ${params.cores} \
      ${genome_file} UMI_collapsed.fq.gz ${bed_arg} | \
    ${params.samtools} view -bS -@ ${params.cores} - | \
    ${params.samtools} sort - -@ ${params.cores} -o UMI_collapsed.bam && \
    ${params.samtools} index UMI_collapsed.bam
    """
}


workflow {
    // Stage 1: Annotation
    if(params.gtf_path || params.gene_bed_path){
        gtf = params.gtf_path ? file(params.gtf_path) : file('NO_FILE')
        gene_bed = params.gene_bed_path ? file(params.gene_bed_path) : file('NO_FILE')
        
        annotation_result = RUN_ANNOTATION(
            gtf,
            gene_bed
        )
    }
    
    // Stage 2: Barcode extraction and mapping
    if(params.fastq_path && params.barcode_path){
        fastq = file(params.fastq_path)
        barcode = file(params.barcode_path)
        
        fastq_chunks = SPLIT_FASTQ(fastq).flatten()
        
        chunk_results = EXTRACT_TAG_BC_CHUNK(
            fastq_chunks,
            barcode
        )
        
        merged_barcode = MERGE_BARCODE_MATCH(
            chunk_results[1].collect()
        )
        
        // Stage 3: Map polished FASTQ chunks to genome (parallel)
        if(params.genome_path){
            genome = file(params.genome_path)
            bed = params.minimap_bed_path ? file(params.minimap_bed_path) : file('NO_FILE')
            
            bam_chunks = MAP_POLISHED_FASTQ_CHUNK(
                chunk_results[0],
                genome,
                bed
            )

            mapping_result = MERGE_BAM(
                bam_chunks.collect()
            )

            // Stage 4: Extract isoforms and integrate barcodes
            if(params.genome_name && (params.gtf_path || params.gene_bed_path)){
                // Filter genes without BAM coverage
                filtered_gene_bed = FILTER_GENE_BED(
                    mapping_result[0],
                    mapping_result[1],
                    annotation_result[0]
                )

                barcode_file_ch = merged_barcode.first()
                gene_bed_ch = filtered_gene_bed.first()

                // Use per-chunk BAMs directly for parallel extract_and_integrate
                chunk_iso_results = EXTRACT_AND_INTEGRATE_READS_CHUNK(
                    barcode_file_ch,
                    bam_chunks,
                    gene_bed_ch
                )

                // 4d: Merge chunk results
                extract_and_integrate_result = MERGE_ISOFORM_READS(
                    chunk_iso_results.collect()
                )

                // Stage 5: UMI deduplication
                // Use .first() so the single file is broadcast to all parallel chunks
                gene_bed_rds = annotation_result[0].first()
                
                split_results = SPLIT_UMI_DATA(
                    extract_and_integrate_result[0],  // BarcodeMatchIso.txt
                    extract_and_integrate_result[1]   // adapterNeedle.txt
                )
                
                chunk_data = split_results[0].flatten()
                chunk_qual = split_results[1].flatten()
                
                chunk_umi_results = UMI_COUNT_PARALLEL_CHUNK(
                    chunk_data,
                    chunk_qual,
                    gene_bed_rds
                )
                
                umi_count_result = MERGE_UMI_COUNT(
                    chunk_umi_results.flatten().collect()
                )

                // Stage 6: Isoform imputation (conditional on to_isoform flag)
                if(params.to_isoform && params.gtf_path){
                    // Use .first() so the single GTF file is broadcast to all parallel chunks
                    gtf_rds = annotation_result[1].first()

                    iso_chunks = SPLIT_ISO_INPUT(umi_count_result).flatten()

                    iso_chunk_results = UMI_COUNT_TO_ISOFORM_CHUNK(
                        iso_chunks,
                        gtf_rds
                    )

                    iso_mat_result = MERGE_ISO_MAT(
                        iso_chunk_results.collect()
                    )
                }

                // Stage 7: UMI consensus output and remapping (chunked)
                consensus_chunks = SPLIT_CONSENSUS_INPUT(umi_count_result).flatten()

                consensus_fastq_results = UMI_CONSENSUS_FASTQ_CHUNK(
                    consensus_chunks,
                    gene_bed_rds
                )

                consensus_result = MERGE_CONSENSUS_AND_MAP(
                    consensus_fastq_results[0].collect(),
                    consensus_fastq_results[1].collect(),
                    genome,
                    bed
                )
            }
        }
    }
}
