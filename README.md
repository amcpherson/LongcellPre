# LongcellPre
LongcellPre is an R pipeline to analyze Nanopore long read sequencing dataset based on 10X single cell sequencing toolkit. This pipeline includes preprocessing to do barcode and unique molecular identifier (UMI) assignment to give an accurate isoform quantification. Based on the isoform quantification from LongcellPre, our another  pipeline Longcell incorporates downstream splicing analysis, including identification of highly variable exons and differential alternative splicing analysis between different cell populations.

## Installation
requires:  
- minimap2: https://github.com/lh3/minimap2
- samtools: http://www.htslib.org/
- bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html

You can install `LongcellPre` by `devtools`:
```
devtools::install_github("yuntianf/Longcellsrc")
devtools::install_github("yuntianf/LongcellPre",dependencies=TRUE)
```
The singularity image for this package is under construction.

## Workflow
The simplist way to run `LongcellPre` is to apply:
```
Rscript ./LongcellPre/exec/RunLongcellPre.R -f {FASTQ} -b {BARCODES} -t {TOOLKIT} -q {PROTOCOL} -g {GENOME} -n {GENOME_NAME} --gtf {GTF} -o {OUTDIR}
```
The execution file can also be downloaded from https://github.com/yuntianf/LongcellPre/blob/main/exec/RunLongcellPre.R


The basic parameters for this pipeline are shown here:
```
Pipeline for single cell Nanopore RNA-seq preprocessing

options:
  -h, --help            show this help message and exit
  -f FASTQ, --fastq FASTQ
                        The path for the input fastq file
  -b BARCODE, --barcode BARCODE
                        The path for the cell barcode whitelist
  -t {5,3}, --toolkit {5,3}
                        The toolkit used in sequencing, should be 5 or 3
  -q {10X,VISIUM,Curio,other}, --protocol {10X,VISIUM,Curio,other}
                        The sequencing protocol, ex. '10X', 'VISIUM', 'Curio'
  -g GENOME_PATH, --genome_path GENOME_PATH
                        The path of the genome reference
  -n GENOME_NAME, --genome_name GENOME_NAME
                        the genome name used for mapping, ex. hg38
  --gtf GTF             The path of the gtf annotation
  --gene_bed_path GENE_BED_PATH
                        The path of the gene bed annotation
  --minimap_bed_path MINIMAP_BED_PATH
                        The path of your bed annotation for minimap2, can be
                        generated from GTF/GFF3 with ‘paftools.js gff2bed
                        anno.gtf’
  -o WORK_DIR, --work_dir WORK_DIR
                        The output directory
  --to_isoform TO_ISOFORM
                        A flag to indicate if the cell by isoform matrix
                        should be generated
  -c CORES, --cores CORES
                        The number of cores used for parallization
  -m {sequential,multisession,multicore,cluster}, --mode {sequential,multisession,multicore,cluster}
                        The mode for parallization. The parallization is
                        implemented with future.apply, the feasible modes can
                        be 'sequential','multicore','cluster'
  --minimap2 MINIMAP2   The path of the minimap2
  --samtools SAMTOOLS   The path of the samtools
  --bedtools BEDTOOLS   The path of the bedtools
```
To view all parameters for this pipeline, you can run `Rscript ./LongcellPre/exec/RunLongcellPre.R -h --full`

We provide a demo data with 244 cells and 3 genes. This data is a subset of the colorectal metastasis sample we used in the paper. The data and corresponding annotations can be downloaded from: 
https://www.dropbox.com/scl/fo/21tw8rrkaancani0fzq3t/AKNHUk06onR2c2dYuB4wXWY?rlkey=1zikug28qr9ziw2cdsgelrm9p&st=ypm9m00i&dl=0

The demo data can be processed by:
```
Rscript RunLongcellPre.R -f example.fastq.gz -b barcodes.txt -t 5 -q 10X -g genome.fa -n hg38 --gtf gencode.v39.sub.gtf -o ./demo/ -c 4 -m multicore
```

## output
The output of the LongcellPre pipeline includes:
```
├── annotation
│   ├── exon_gtf.rds: The exon annotation for each gene given the gtf annotation, stored in RDS format.
│   ├── exon_gtf.txt: The exon annotation for each gene given the gtf annotation, stored in tsv format.
│   ├── gene_bed.rds: The annotation non-overlapped exon bins for each gene given the gtf annotation, stored in RDS format.
│   ├── gene_bed.txt: The annotation non-overlapped exon bins for each gene given the gtf annotation, stored in tsv format.
├── bam
│   ├── polish.bam: The mapped bam from reads with adapters trimmed.
│   └── polish.bam.bai
├── BarcodeMatch
│   ├── adapterNeedle.txt: The summary statistics for the Needleman score of the kmer adapter sequence in each reads compared to its know sequence.
│   └── BarcodeMatchIso.txt: The cell barcode, UMI, and read alignment for each read.
├── out
│   ├── gene: The sparse matrix for the cell by gene matrix.
│   │   ├── barcodes.tsv
│   │   ├── features.tsv
│   │   └── matrix.mtx
│   ├── iso_count.txt: The isoform quantification in each cell (without isoform annotation).
│   ├── isoform: The sparse matrix for the cell by isoform matrix.
│   │   ├── barcodes.tsv
│   │   ├── features.tsv
│   │   └── matrix.mtx
│   ├── reads_annot.csv: The read annotation for the UMI collapsed fastq file.
│   └── UMI_collapsed.fq.gz: The polished reads after UMI collapsion. This can be input for other tools.
│   └── UMI_collapsed.bam: The mapping result for UMI_collapsed.fq.gz.
│   └── UMI_collapsed.bam.bai
├── arg.log: The log file to record used parameters.
└── polish.fq.gz
```
The isoform quantification output `iso_count.txt` can be directly analyzed by our downstream package `Longcell` (https://github.com/yuntianf/Longcell).

For the isoform quantification formatted as cell by isoform matrix, we recommend the combination of LongcellPre and IsoQuant (https://github.com/ablab/IsoQuant), which has the overall best performance in our benchmark. The input for IsoQuant is `UMI_collapsed.bam` and `reads_annot.csv`:
```
python isoquant.py --reference $GENOME_PATH \
--genedb $GTF \
--bam UMI_collapsed.bam \
--read_group file:reads_annot.csv:0:1:, \
--data_type nanopore -o $OUTDIR \
--report_novel_unspliced true \
--clean_start
```

## Nextflow Usage

The LongcellPre pipeline is also available as a Nextflow DSL2 workflow for scalable, containerized execution. This version orchestrates all 7 stages of the pipeline with automatic dependency resolution and built-in parallelization.

### Requirements

- **Nextflow**: version 24.04.4 or later
- **Docker** or **Singularity**: for containerized execution
- **Optional**: micromamba for environment management

### Quick Start

Run directly from GitHub without cloning:

```bash
# Basic command with GTF annotation
nextflow run yuntianf/LongcellPre \
  --gtf_path <path_to_gtf> \
  --fastq_path <path_to_fastq> \
  --barcode_path <path_to_barcodes> \
  --genome_path <path_to_genome> \
  --genome_name <BSgenome_name> \
  --cores 4

# With gene BED instead of GTF
nextflow run yuntianf/LongcellPre \
  --gene_bed_path <path_to_gene_bed> \
  --fastq_path <path_to_fastq> \
  --barcode_path <path_to_barcodes> \
  --genome_path <path_to_genome> \
  --genome_name <BSgenome_name> \
  --results_dir ./output \
  --cores 4
```

### Nextflow Parameters

Key parameters for the Nextflow workflow:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gtf_path` | null | Path to GTF annotation file |
| `gene_bed_path` | null | Path to gene BED annotation (alternative to GTF) |
| `fastq_path` | null | **Required**: Path to input FASTQ file |
| `barcode_path` | null | **Required**: Path to cell barcode whitelist |
| `genome_path` | null | **Required**: Path to reference genome FASTA |
| `genome_name` | null | **Required**: BSgenome name (e.g., `BSgenome.Hsapiens.UCSC.hg38`) |
| `minimap_bed_path` | null | BED file with splicing sites for minimap2 |
| `toolkit` | 5 | Sequencing toolkit: 5 or 3 |
| `protocol` | '10X' | Sequencing protocol (10X, VISIUM, Curio, other) |
| `cores` | 4 | Number of CPU cores for parallelization |
| `results_dir` | 'results' | Output directory for results |
| `to_isoform` | true | Generate cell-by-isoform matrix (requires GTF) |
| `map_qual` | 30 | Mapping quality threshold |
| `end_flank` | 200 | End flanking region for reads |
| `splice_site_bin` | 2 | Splice site bin size |
| `fastq_chunk_size` | 1000000 | Number of reads per FASTQ chunk for barcode extraction |
| `genes_per_chunk` | 500 | Approximate number of genes per chunk for UMI counting and isoform imputation |

### Example with Test Data

```bash
# Run with demo data (from test_data/)
nextflow run yuntianf/LongcellPre \
  --gtf_path test_data/gencode.v39.sub.gtf \
  --fastq_path test_data/example.fq.gz \
  --barcode_path test_data/barcodes.txt \
  --genome_path test_data/genome.fa \
  --genome_name BSgenome.Hsapiens.UCSC.hg38 \
  --results_dir ./demo_output \
  --cores 4
```

### Pipeline Stages

The Nextflow workflow implements all 7 stages of LongcellPre, with built-in chunking for parallelization:

1. **RUN_ANNOTATION** - Generate gene annotation from GTF/BED
2. **EXTRACT_TAG_BC** - Split FASTQ into chunks (`fastq_chunk_size` reads each), extract barcodes/UMI from each chunk in parallel, then merge polished FASTQ and barcode results
3. **MAP_POLISHED_FASTQ** - Map polished reads to genome with minimap2
4. **EXTRACT_AND_INTEGRATE_READS** - Filter genes with BAM coverage, split gene BED into chunks (`genes_per_chunk`), extract isoforms and integrate barcodes for each chunk in parallel, then merge
5. **UMI_COUNT** - Split barcode/isoform data into gene chunks, run UMI deduplication for each chunk in parallel, then merge
6. **UMI_COUNT_TO_ISOFORM** - Split UMI counts into gene chunks, run isoform imputation for each chunk in parallel, then merge (optional, requires GTF)
7. **UMI_CONSENSUS_OUT** - Generate UMI-collapsed FASTQ and BAM

All stages run in Docker containers for reproducibility and are automatically orchestrated by Nextflow. Chunking parameters (`fastq_chunk_size`, `genes_per_chunk`) control the granularity of parallelization.

## Citation

If you use Longcell for published work, please cite our manuscript:

``` r
Fu Y, Kim H, Roy S, et al. Single cell and spatial alternative splicing analysis with Nanopore long read sequencing[J]. Nature Communications, 2025, 16(1): 6654.
```
