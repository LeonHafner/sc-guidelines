FROM bioconductor/bioconductor_docker:RELEASE_3_16

# Install the remotes package
RUN R -e "install.packages('remotes')"

# Install older version of ggplot2 and cowplot to circumvent this issue: https://github.com/wilkelab/cowplot/issues/202
# (use dependencies=F to prevent overwrite of ggplot2 version to 3.5)
RUN R -e "remotes::install_version('ggplot2', version='3.4.4', repos='http://cran.rstudio.com/', dependencies = T)"
RUN R -e "remotes::install_version('cowplot', version='1.1.2', repos='http://cran.rstudio.com/', dependencies = F)"

RUN R -e "install.packages(c('argparser', 'data.table', 'ggokabeito', 'pracma', 'dplyr', 'png'))"

# Install anndataR and dependencies
RUN R -e "BiocManager::install(c('rhdf5', 'SingleCellExperiment'))"
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('scverse/anndataR')"
