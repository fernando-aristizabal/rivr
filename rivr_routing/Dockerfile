FROM rocker/r-devel:latest

# INSTALL PACKAGES
RUN R -e "install.packages('rivr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('doParallel',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('rhdf5')"

# USER AND GROUP
ARG GroupID=1370800235
ARG GroupName=fim
RUN addgroup --gid $GroupID $GroupName
USER root:$GroupName
