FROM bioconductor/bioconductor_docker:RELEASE_3_16

RUN R -e "BiocManager::install(c('distinct', 'rhdf5', 'SingleCellExperiment', 'scuttle'))"
RUN R -e "install.packages(c('argparser', 'devtools'))"

RUN R -e "devtools::install_github('scverse/anndataR')"