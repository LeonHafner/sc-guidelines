FROM bioconductor/bioconductor_docker:RELEASE_3_16

RUN R -e "BiocManager::install(c('splatter', 'scater', 'VariantAnnotation', 'rhdf5', 'SingleCellExperiment', 'MAST'))"
RUN R -e "install.packages(c('argparser', 'devtools', 'lme4'))"

RUN git clone https://github.com/bmbolstad/preprocessCore.git
RUN R CMD INSTALL --configure-args="--disable-threading" preprocessCore

RUN R -e "devtools::install_github('scverse/anndataR')"

