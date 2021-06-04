FROM rocker/r-devel:latest

RUN R -e "install.packages('rivr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('doParallel',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('rhdf5')"

#RUN apt update && \
#    apt install -y postgresql \
#    postgresql-contrib \
#    && rm -rf /var/lib/apt/lists/*

