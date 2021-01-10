FROM rocker/rstudio:3.6.3

RUN install2.r --error \
    --deps TRUE \
    devtools \
    rlang \
    optparse \
    docstring \
    openxlsx \
    Rserve \
    RColorBrewer \
    xtable \
    som \
    ROCR \
    RJSONIO \
    gplots \
    e1071 \
    caTools \
    igraph \
    randomForest \
    Cairo \
    pls \
    pheatmap \
    lattice \
    rmarkdown \
    knitr \
    data.table \
    pROC \
    Rcpp \
    caret \
    ellipse
    
RUN R -e "BiocManager::install('xcms')"
RUN R -e "BiocManager::install('impute')"
RUN R -e "BiocManager::install('pcaMethods')"
RUN R -e "BiocManager::install('siggenes')"
RUN R -e "BiocManager::install('globaltest')"
RUN R -e "BiocManager::install('GlobalAncova')"
RUN R -e "BiocManager::install('Rgraphviz')"
RUN R -e "BiocManager::install('pcaMethods')"
RUN R -e "BiocManager::install('siggenes')"
RUN R -e "BiocManager::install('KEGGgraph')"
RUN R -e "BiocManager::install('preprocessCore')"
RUN R -e "BiocManager::install('genefilter')"
RUN R -e "BiocManager::install('SSPA')"
RUN R -e "BiocManager::install('sva')"

RUN R -e "devtools::install_github('xia-lab/MetaboAnalystR', build = TRUE, build_vignettes = TRUE, build_manual = TRUE)"

ADD ./ /scripts
ADD ./ /notebooks
ADD ./ /data_examples
