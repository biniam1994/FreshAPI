FROM rocker/r-base

RUN apt-get update -qq && apt-get install -y \
  bash \
  curl \
  pandoc \
  pandoc-citeproc \
  git-core \
  libssl-dev \
  libcurl4-gnutls-dev


RUN R -e "install.packages(c('plumber', 'caret', 'randomForest', 'ggplot2', 'ggpubr', 'zoo', 'jsonlite', 'stats', 'kernlab', 'rmarkdown', 'tinytex', 'openxlsx', 'prospectr', 'readxl', 'factoextra', 'dplyr', 'gmodels', 'assertthat'))"

RUN R -e "tinytex::install_tinytex()"

COPY plumber.R /plumber.R
COPY plumber-api_router.R /plumber-api_router.R
COPY functions.R /functions.R
COPY models /models
COPY testData /testData

WORKDIR /freshAPI_GP

ENTRYPOINT ["R", "-e", "source('plumber-api_router.R')"]




