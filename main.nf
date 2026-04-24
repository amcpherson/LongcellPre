nextflow.enable.dsl=2

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
    path 'annotation'

    script:
    def gtf_arg = gtf_file.name != 'gene_bed_file' ? "--gtf_path ${gtf_file}" : ''
    def gene_bed_arg = gene_bed_file.name != 'gene_bed_file' ? "--gene_bed_path ${gene_bed_file}" : ''
    def overwrite_arg = params.overwrite ? 'TRUE' : 'FALSE'
    """
    mkdir -p annotation
    annotation.R ${gtf_arg} ${gene_bed_arg} --work_dir . --overwrite ${overwrite_arg}
    """
}

process EXTRACT_TAG_BC {
    publishDir "${params.work_dir}", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path fastq_file
    path barcode_file

    output:
    path 'polish.fq.gz'
    path 'BarcodeMatch.txt'

    script:
    def adapter_arg = params.adapter ? "--adapter ${params.adapter}" : ''
    """
    extract_tag_bc.R \
      --fastq_path ${fastq_file} \
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
    path annotation_dir

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
      --annotation_dir ${annotation_dir} \
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

process UMI_COUNT_PARALLEL {
    publishDir "${params.work_dir}/out", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path barcode_iso_file
    path adapter_needle_file
    path gene_bed_file

    output:
    path 'iso_count.txt'

    script:
    def sim_thresh_arg = params.sim_thresh ? "--sim_thresh ${params.sim_thresh}" : ''
    def verbose_arg = params.verbose ? "--verbose TRUE" : "--verbose FALSE"
    """
    mkdir -p out
    umi_count_parallel.R \
      --data_path ${barcode_iso_file} \
      --qual_path ${adapter_needle_file} \
      --gene_bed_path ${gene_bed_file} \
      --out_dir . \
      --splice_site_thresh ${params.splice_site_thresh} \
      ${sim_thresh_arg} \
      ${verbose_arg} \
      --bed_gene_col ${params.bed_gene_col} \
      --bed_strand_col ${params.bed_strand_col} \
      --cores ${params.cores}
    """
}

process UMI_COUNT_TO_ISOFORM {
    publishDir "${params.work_dir}/out", mode: 'copy', overwrite: true
    container 'quay.io/andrew_mcpherson/longcellpre:latest'

    input:
    path umi_count_file
    path gene_bed_file
    path gtf_file

    output:
    path 'iso_count_mat.txt', optional: true

    script:
    def gtf_arg = gtf_file.name != 'NO_FILE' ? "--gtf_path ${gtf_file}" : "--gtf_path NO_FILE"
    def filter_arg = params.filter_only_intron ? "--filter_only_intron TRUE" : "--filter_only_intron FALSE"
    """
    mkdir -p out
    umi_count_to_isoform.R \
      --umi_count_path ${umi_count_file} \
      --gene_bed_path ${gene_bed_file} \
      ${gtf_arg} \
      --out_dir . \
      ${filter_arg} \
      --mid_offset_thresh ${params.mid_offset_thresh} \
      --overlap_thresh ${params.overlap_thresh} \
      --gtf_gene_col ${params.gtf_gene_col} \
      --gtf_start_col ${params.gtf_start_col} \
      --gtf_end_col ${params.gtf_end_col} \
      --gtf_iso_col ${params.gtf_iso_col} \
      --split ${params.split} \
      --sep ${params.sep} \
      --bed_gene_col ${params.bed_gene_col} \
      --bed_strand_col ${params.bed_strand_col} \
      --cores ${params.cores}
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
        
        extract_result = EXTRACT_TAG_BC(
            fastq,
            barcode
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
                    annotation_result
                )

                // Stage 5: UMI deduplication
                gene_bed_rds = file("${params.work_dir}/annotation/gene_bed.rds")
                umi_count_result = UMI_COUNT_PARALLEL(
                    extract_and_integrate_result[0],  // BarcodeMatchIso.txt
                    extract_and_integrate_result[1],  // adapterNeedle.txt
                    gene_bed_rds
                )

                // Stage 6: Isoform imputation (conditional on to_isoform flag)
                if(params.to_isoform && params.gtf_path){
                    gtf_rds = file("${params.work_dir}/annotation/exon_gtf.rds")
                    iso_impute_result = UMI_COUNT_TO_ISOFORM(
                        umi_count_result,
                        gene_bed_rds,
                        gtf_rds
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
