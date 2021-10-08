FROM rocker/shiny:4.1.0
LABEL authors="Alex Lemenze" \
    description="Docker image containing the scRNA template."

RUN apt-get update && apt-get install -y \ 
    sudo libhdf5-dev build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

RUN R -e "install.packages(c('Seurat','DT','RColorBrewer','colourpicker','tidyverse','shinyjs','shinythemes','shinycssloaders','hdf5r','dplyr','cowplot','knitr','slingshot','msigdbr','remotes','metap','devtools','R.utils','stringr'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('SingleR','slingshot','BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'Matrix.utils','scRNAseq','celldex','fgsea','multtest','scuttle'))"
RUN R -e "remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')"
RUN R -e "remotes::install_github('satijalab/seurat-wrappers')"
RUN R -e "devtools::install_github('cole-trapnell-lab/monocle3')"
RUN R -e "BiocManager::install(c('clusterProfiler','AnnotationHub','enrichplot'))"
RUN R -e "install.packages('rhandsontable','circlize')"

COPY ./app /srv/shiny-server/
COPY shiny-customized.config /etc/shiny-server/shiny-server.conf
RUN sudo chown -R shiny:shiny /srv/shiny-server
EXPOSE 8080

USER shiny
CMD ["/usr/bin/shiny-server"]