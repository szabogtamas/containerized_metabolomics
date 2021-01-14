FROM rocker/rstudio:3.6.3

RUN sudo apt-get update -y
RUN sudo apt-get install -y libxt-dev
RUN sudo apt-get install -y libx11-dev
RUN sudo apt-get install -y libmagick++-dev
RUN sudo apt-get install -y libfontconfig1-dev
RUN sudo apt-get install -y zlib1g-dev
RUN sudo apt-get install -y cargo
RUN sudo apt-get install -y libxml2-dev
RUN sudo apt-get install -y libcairo2-dev
RUN sudo apt-get install -y mesa-common-dev
RUN sudo apt-get install -y libglu1-mesa-dev
RUN sudo apt-get install -y libopenmpi-dev

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
