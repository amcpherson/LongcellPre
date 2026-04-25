# Use rocker/r-ver as base image with R 4.4.0
FROM rocker/r-ver:4.4.0

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies commonly needed for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libglpk-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools
RUN apt-get update && apt-get install -y \
    minimap2 \
    samtools \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# Install remotes package (lighter than devtools, sufficient for GitHub installs)
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org/', Ncpus=2)" \
    && R -e "if (!require('remotes')) stop('remotes installation failed')"

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/', Ncpus=2)" \
    && R -e "if (!require('BiocManager')) stop('BiocManager installation failed')"

# Install common CRAN dependencies
RUN R -e "install.packages(c('matrixStats', 'Matrix', 'Rcpp', 'RcppArmadillo'), repos='https://cloud.r-project.org/', Ncpus=2)"

# Pre-install critical Bioconductor packages in dependency order
RUN R -e "BiocManager::install(c('BiocGenerics', 'S4Vectors', 'IRanges', 'GenomeInfoDb'), ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install(c('Rhtslib', 'Rsamtools', 'zlibbioc'), ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install(c('XVector', 'Biostrings', 'GenomicRanges', 'GenomicAlignments'), ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install(c('rtracklayer', 'ShortRead', 'BSgenome'), ask=FALSE, update=FALSE)"

# Install LongcellPre dependencies
RUN R -e "remotes::install_github('yuntianf/Longcellsrc')"
RUN R -e "remotes::install_github('yuntianf/LongcellPre', dependencies=TRUE)"

# Clone the LongcellPre repository to access exec/ scripts
RUN git clone https://github.com/yuntianf/LongcellPre.git /opt/LongcellPre

# Add exec scripts to PATH
ENV PATH="/opt/LongcellPre/exec:${PATH}"

# Set working directory
WORKDIR /workspace

# Default command
CMD ["R"]
