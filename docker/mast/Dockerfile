FROM bioconductor/bioconductor_docker:RELEASE_3_16

RUN R -e "BiocManager::install(c('MAST', 'rhdf5', 'SingleCellExperiment'))"
RUN R -e "install.packages(c('argparser', 'devtools', 'lme4'))"

RUN R -e "devtools::install_github('scverse/anndataR')"