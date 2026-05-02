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
// The barcode/isoform data is split into chunks based on genes (~50 genes per
// chunk, configurable via params.genes_per_chunk). Each gene chunk is processed
// independently for UMI deduplication, then results are merged.
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
params.work_dir = 'test_data/nextflow_annotation_output'
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
params.mu = 20
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
params.genes_per_chunk = 50  // Approximate number of genes per chunk
// Tool paths
params.minimap2 = "minimap2"
params.samtools = "samtools"

process RUN_ANNOTATION {
    publishDir "${params.work_dir}", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path gtf_file
    path gene_bed_file

    output:
    path 'annotation/gene_bed.rds'
    path 'annotation/exon_gtf.rds', optional: true

    script:
    def gtf_arg = gtf_file.name != 'gene_bed_file' ? "--gtf_path ${gtf_file}" : ''
    def gene_bed_arg = gene_bed_file.name != 'gene_bed_file' ? "--gene_bed_path ${gene_bed_file}" : ''
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
    gunzip -c '${fastq_file}' | split -l \$((${params.fastq_chunk_size} * 4)) --numeric-suffixes=1 -d - chunk_raw_
    
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

process MERGE_POLISH_FASTQ {
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path 'polish_*.fq.gz'
    path 'barcode_*.txt'

    output:
    path 'polish.fq.gz'
    path 'BarcodeMatch.txt'

    publishDir "${params.work_dir}", mode: 'copy', overwrite: true

    script:
    """
    # Merge polish FASTQ files
    cat polish_*.fq.gz > polish.fq.gz
    
    # Merge barcode match files, keeping header from first file
    awk 'FNR==1 && NR>1 { next; } { print }' barcode_*.txt > BarcodeMatch.txt
    """
}

process MAP_POLISHED_FASTQ {
    publishDir "${params.work_dir}/bam", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path fastq_file
    path genome_file
    path bed_file

    output:
    path 'polish.bam'
    path 'polish.bam.bai'

    script:
    def bed_arg = bed_file.name != 'NO_FILE' ? "--junc-bed ${bed_file}" : ''
    """
    minimap2 -ax splice -uf --sam-hit-only -t ${params.cores} ${bed_arg} ${genome_file} ${fastq_file} | \
    samtools view -bS -@ ${params.cores} - | \
    samtools sort - -@ ${params.cores} -o polish.bam && \
    samtools index polish.bam
    """
}

process EXTRACT_AND_INTEGRATE_READS {
    publishDir "${params.work_dir}", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path fastq_file
    path barcode_file
    path bam_file
    path bai_file
    path gene_bed_rds

    output:
    path 'BarcodeMatch/BarcodeMatchIso.txt'
    path 'BarcodeMatch/adapterNeedle.txt'

    script:
    def adapter_arg = params.adapter ? "--adapter ${params.adapter}" : ''
    def bed_arg = params.minimap_bed_path ? "--minimap_bed_path ${params.minimap_bed_path}" : ''
    """
    mkdir -p BarcodeMatch
    cp ${barcode_file} BarcodeMatch/BarcodeMatch.txt
    mkdir -p bam
    cp ${bam_file} bam/polish.bam
    cp ${bai_file} bam/polish.bam.bai
    extract_and_integrate_reads.R \
      --fastq_path ${fastq_file} \
      --barcode_path BarcodeMatch/BarcodeMatch.txt \
      --gene_bed_path ${gene_bed_rds} \
      --work_dir . \
      --genome_path ${params.genome_path} \
      --genome_name ${params.genome_name} \
      --toolkit ${params.toolkit} \
      --protocol ${params.protocol} \
      ${adapter_arg} \
      ${bed_arg} \
      --minimap2 ${params.minimap2} \
      --samtools ${params.samtools} \
      --bedtools ${params.bedtools} \
      --map_qual ${params.map_qual} \
      --end_flank ${params.end_flank} \
      --splice_site_bin ${params.splice_site_bin} \
      --mean_edit_thresh ${params.mean_edit_thresh} \
      --cores ${params.cores}
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
    data <- read.table('${barcode_iso_file}', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    qual <- read.table('${adapter_needle_file}', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Get unique genes and split into chunks
    genes <- unique(data\$gene)
    n_genes <- length(genes)
    n_chunks <- max(1, ceiling(n_genes / ${params.genes_per_chunk}))
    
    cat("Splitting", n_genes, "genes into", n_chunks, "chunks\\n")
    
    # Create gene chunks - roughly balanced
    gene_chunks <- split(genes, ceiling(seq_along(genes) / (n_genes / n_chunks)))
    
    # Write each chunk
    for (i in seq_along(gene_chunks)) {
        chunk_genes <- gene_chunks[[i]]
        chunk_data <- data[data\$gene %in% chunk_genes, ]
        
        chunk_file <- paste0('data_chunk_', sprintf('%03d', i), '.txt')
        write.table(chunk_data, file = chunk_file, sep = "\t", quote = FALSE, 
                    row.names = FALSE, col.names = TRUE)
        cat("Wrote", nrow(chunk_data), "rows to", chunk_file, "\\n")
    }
    
    # Qual file is the same for all chunks (no need to split), but copy for each chunk
    for (i in seq_along(gene_chunks)) {
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
    publishDir "${params.work_dir}/out", mode: 'copy', overwrite: true
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

    data <- read.table('${umi_count_file}', header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    genes <- unique(data\$gene)
    n_genes <- length(genes)
    n_chunks <- max(1, ceiling(n_genes / ${params.genes_per_chunk}))

    cat("Splitting", n_genes, "genes into", n_chunks, "chunks\\n")

    gene_chunks <- split(genes, ceiling(seq_along(genes) / (n_genes / n_chunks)))

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
    publishDir "${params.work_dir}/out", mode: 'copy', overwrite: true
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
        read.table(f, header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
    })
    combined <- as.data.frame(do.call(rbind, chunks))
    # Drop empty-header rows from chunks with no results
    combined <- combined[nchar(combined\$cell) > 0, ]
    cat("Combined", nrow(combined), "rows\\n")

    saveIsoMat(combined, ".")
    cat("Isoform matrix saved\\n")
    """
}

process UMI_CONSENSUS_OUT {
    publishDir "${params.work_dir}/out", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path umi_count_file
    path gene_bed_file
    path genome_file
    path bed_file

    output:
    path 'UMI_collapsed.fq.gz'
    path 'reads_annot.csv'
    path 'UMI_collapsed.bam'

    script:
    def bed_arg = bed_file.name != 'NO_FILE' ? "--minimap_bed_path ${bed_file}" : ''
    """
    mkdir -p out
    umi_consensus_out.R \
      --umi_count_path ${umi_count_file} \
      --gene_bed_path ${gene_bed_file} \
      --genome_path ${genome_file} \
      --genome_name ${params.genome_name} \
      ${bed_arg} \
      --minimap2 ${params.minimap2} \
      --samtools ${params.samtools} \
      --out_dir . \
      --cores ${params.cores}
    """
}


workflow {
    // Stage 1: Annotation
    if(params.gtf_path || params.gene_bed_path){
        gtf = file(params.gtf_path)
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
        
        extract_result = MERGE_POLISH_FASTQ(
            chunk_results[0].collect(),
            chunk_results[1].collect()
        )
        
        // Stage 3: Map polished FASTQ to genome
        if(params.genome_path){
            genome = file(params.genome_path)
            bed = params.minimap_bed_path ? file(params.minimap_bed_path) : file('NO_FILE')
            
            mapping_result = MAP_POLISHED_FASTQ(
                extract_result[0],  // polish.fq.gz
                genome,
                bed
            )

            // Stage 4: Extract isoforms and integrate barcodes
            if(params.genome_name && (params.gtf_path || params.gene_bed_path)){
                extract_and_integrate_result = EXTRACT_AND_INTEGRATE_READS(
                    extract_result[0],
                    extract_result[1],
                    mapping_result[0],
                    mapping_result[1],
                    annotation_result[0]
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

                // Stage 7: UMI consensus output and remapping
                consensus_result = UMI_CONSENSUS_OUT(
                    umi_count_result,
                    gene_bed_rds,
                    genome,
                    bed
                )
            }
        }
    }
}
