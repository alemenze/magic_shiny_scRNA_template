FROM rocker/shiny:4.1.0
LABEL authors="Alex Lemenze" \
    description="Docker image containing the scRNA template."

RUN apt-get update && apt-get install -y \ 
    sudo libssl-dev libv8-dev libsodium-dev libxml2-dev \
    libcurl4-openssl-dev libhdf5-dev

RUN R -e "install.packages(c('shiny','DT','jpeg','readr','ggthemes','markdown','tidyr','rmarkdown','Hmisc','tidyverse','r-lib','shinythemes','shinycssloaders','RColorBrewer','pheatmap','ggplot2','shinyjs','dplyr','colourpicker'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('GenomicRanges','SummarizedExperiment','Biobase','genefilter','locfit','geneplotter','Hmisc','RcppArmadillo','msigdbr','stringr'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('remotes'),repos='http://cran.rstudio.com/')"
RUN R -e 'install.packages("XML", repos = "http://www.omegahat.net/R")'
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('annotate','geneplotter','genefilter','clusterProfiler','pathview','enrichplot'))"
RUN R -e "BiocManager::install('org.Mm.eg.db', character.only = TRUE)"
RUN R -e "BiocManager::install('org.Hs.eg.db', character.only = TRUE)"
RUN R -e "install.packages(c('httpuv'),repos='http://cran.rstudio.com/')"

RUN R -e "install.packages(c('Seurat','hdf5r','dplyr','cowplot','knitr','slingshot','msigdbr','metap','devtools','R.utils'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('SingleR','slingshot','scRNAseq','celldex','fgsea','multtest','scuttle','BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'Matrix.utils'))"
RUN R -e "remotes::install_github('satijalab/seurat-wrappers')"
RUN R -e "devtools::install_github('cole-trapnell-lab/monocle3')"

COPY ./app /srv/shiny-server/
COPY shiny-customized.config /etc/shiny-server/shiny-server.conf
RUN sudo chown -R shiny:shiny /srv/shiny-server
EXPOSE 8080

USER shiny
CMD ["/usr/bin/shiny-server"]