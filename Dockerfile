FROM rocker/verse:4.1

RUN sudo apt-get update -y &&\
  sudo apt-get install -y libxt-dev &&\
  sudo apt-get install -y libx11-dev &&\
  sudo apt-get install -y curl &&\
  sudo apt-get install -y libmagick++-dev &&\
  sudo apt-get install -y libfontconfig1-dev &&\
  sudo apt-get install -y zlib1g-dev &&\
  sudo apt-get install -y cargo &&\
  sudo apt-get install -y libxml2-dev &&\
  sudo apt-get install -y libcairo2-dev &&\
  sudo apt-get install -y mesa-common-dev &&\
  sudo apt-get install -y libglu1-mesa-dev &&\
  sudo apt-get install -y libopenmpi-dev &&\
  sudo apt-get install -y libnetcdf-dev &&\
  sudo apt-get install -y libssl-dev &&\
  sudo apt-get install -y libcurl4-openssl-dev &&\
  sudo apt-get install -y libssh2-1-dev &&\
  sudo apt-get install -y libpq-dev &&\
  sudo apt-get install -y zlib1g-dev &&\
  sudo apt-get install -y libgeos-dev &&\
  sudo apt-get install -y libudunits2-dev &&\
  sudo apt-get install -y libgdal-dev &&\
  sudo apt-get install -y libsodium-dev

RUN sudo apt-get install -y build-essential chrpath libxft-dev &&\
  sudo apt-get install -y libfreetype6-dev libfreetype6 libfontconfig1-dev libfontconfig1 &&\
  sudo wget https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2 &&\
  sudo tar xvjf phantomjs-2.1.1-linux-x86_64.tar.bz2 -C /usr/local/share/ &&\
  sudo ln -s /usr/local/share/phantomjs-2.1.1-linux-x86_64/bin/phantomjs /usr/local/bin/

RUN cd /tmp &&\
  wget -qO- https://get.nextflow.io | bash &&\
  mv nextflow /usr/local/bin/nextflow  &&\
  sudo chmod 777 /usr/local/bin/nextflow &&\
  sudo chown rstudio /usr/local/bin/nextflow
ENV PATH="/usr/local/bin:${PATH}"

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
    data.table \
    knitr
    
RUN R -e "BiocManager::install('xcms')" &&\
  R -e "BiocManager::install('impute')" &&\
  R -e "BiocManager::install('pcaMethods')" &&\
  R -e "BiocManager::install('siggenes')" &&\
  R -e "BiocManager::install('globaltest')" &&\
  R -e "BiocManager::install('GlobalAncova')" &&\
  R -e "BiocManager::install('Rgraphviz')" &&\
  R -e "BiocManager::install('KEGGgraph')" &&\
  R -e "BiocManager::install('preprocessCore')" &&\
  R -e "BiocManager::install('genefilter')" &&\
  R -e "BiocManager::install('SSPA')" &&\
  R -e "BiocManager::install('BiocParallel')" &&\
  R -e "BiocManager::install('MSnbase')" &&\
  R -e "BiocManager::install('RBGL')" &&\
  R -e "BiocManager::install('multtest')" &&\
  R -e "BiocManager::install('edgeR')" &&\
  R -e "BiocManager::install('limma')" &&\
  R -e "BiocManager::install('fgsea')" &&\
  R -e "BiocManager::install('sva')" &&\
  R -e "BiocManager::install('crmn')" &&\
  R -e "BiocManager::install('ctc')" &&\
  R -e "BiocManager::install('ppcor')" &&\
  R -e "BiocManager::install('graph')" &&\
  R -e "BiocManager::install('vdiffr')"

RUN install2.r --error \
    --deps TRUE \
    igraph \
    randomForest \
    pls \
    pROC \
    Rcpp \
    pryr \
    pkgcond \ 
    glasso \
    huge \
    plotly \
    ggsci \
    readxl \
    googledrive

RUN R -e "devtools::install_github('xia-lab/OptiLCMS')"
RUN R -e "devtools::install_github('xia-lab/MetaboAnalystR', build = TRUE, build_vignettes = TRUE, build_manual = TRUE)"
RUN R -e "tinytex::parse_packages(text = '! LaTeX Error: File `amsmath.sty` not found.')"

RUN sudo add-apt-repository ppa:malteworld/ppa &&\
  sudo apt update  -y &&\
  sudo apt install pdftk

ADD ./ /home/rstudio/repo_files
ADD ./.Rprofile /home/rstudio/.Rprofile
ENV R_PROFILE_USER /home/rstudio/.Rprofile
RUN chmod a+rwx -R /home/rstudio