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
RUN sudo apt-get install -y libnetcdf-dev
RUN sudo apt-get install -y libssl-dev 
RUN sudo apt-get install -y libcurl4-openssl-dev
RUN sudo apt-get install -y libssh2-1-dev
RUN sudo apt-get install -y libpq-dev
RUN sudo apt-get install -y zlib1g-dev
RUN sudo apt-get install -y libgeos-dev
RUN sudo apt-get install -y libudunits2-dev
RUN sudo apt-get install -y libgdal-dev
RUN sudo apt-get install -y libsodium-dev

RUN sudo apt-get install -y build-essential chrpath libxft-dev
RUN sudo apt-get install -y libfreetype6-dev libfreetype6 libfontconfig1-dev libfontconfig1
RUN sudo wget https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2
RUN sudo tar xvjf phantomjs-2.1.1-linux-x86_64.tar.bz2 -C /usr/local/share/
RUN sudo ln -s /usr/local/share/phantomjs-2.1.1-linux-x86_64/bin/phantomjs /usr/local/bin/

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
    Cairo \
    pheatmap \
    lattice \
    rmarkdown \
    knitr \
    data.table
    
RUN R -e "BiocManager::install('xcms')"
RUN R -e "BiocManager::install('impute')"
RUN R -e "BiocManager::install('pcaMethods')"
RUN R -e "BiocManager::install('siggenes')"
RUN R -e "BiocManager::install('globaltest')"
RUN R -e "BiocManager::install('GlobalAncova')"
RUN R -e "BiocManager::install('Rgraphviz')"
RUN R -e "BiocManager::install('KEGGgraph')"
RUN R -e "BiocManager::install('preprocessCore')"
RUN R -e "BiocManager::install('genefilter')"
RUN R -e "BiocManager::install('SSPA')"
RUN R -e "BiocManager::install('BiocParallel')"
RUN R -e "BiocManager::install('MSnbase')"
RUN R -e "BiocManager::install('RBGL')"
RUN R -e "BiocManager::install('multtest')"
RUN R -e "BiocManager::install('edgeR')"
RUN R -e "BiocManager::install('limma')"
RUN R -e "BiocManager::install('fgsea')"
RUN R -e "BiocManager::install('sva')"
RUN R -e "BiocManager::install('crmn')"
RUN R -e "BiocManager::install('ctc')"
RUN R -e "BiocManager::install('ppcor')"
RUN R -e "BiocManager::install('graph')"
RUN R -e "BiocManager::install('vdiffr')"

RUN install2.r --error \
    --deps TRUE \
    igraph \
    randomForest \
    pls \
    pROC \
    Rcpp \
    glasso \
    huge \
    plotly \
    readxl \
    googledrive

RUN R -e "devtools::install_github('xia-lab/MetaboAnalystR', build = TRUE, build_vignettes = TRUE, build_manual = TRUE)"

ADD ./ /home/rstudio
RUN rm -rf /home/rstudio/kitematic
RUN chmod a+rwx -R /home/rstudio
