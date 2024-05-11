FROM bioconductor/bioconductor_docker:RELEASE_3_16

RUN R -e "BiocManager::install(c('variancePartition', 'rhdf5', 'SingleCellExperiment', 'edgeR'))"
RUN R -e "install.packages(c('argparser', 'devtools', 'data.table'))"

RUN R -e "devtools::install_github('scverse/anndataR')"