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
    """
    extract_tag_bc.R --fastq_path ${fastq_file} --barcode_path ${barcode_file}
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
    """
    mkdir -p out
    umi_count_parallel.R \
      --data_path ${barcode_iso_file} \
      --qual_path ${adapter_needle_file} \
      --gene_bed_path ${gene_bed_file} \
      --out_dir . \
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
    """
    mkdir -p out
    umi_count_to_isoform.R \
      --umi_count_path ${umi_count_file} \
      --gene_bed_path ${gene_bed_file} \
      ${gtf_arg} \
      --out_dir . \
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
