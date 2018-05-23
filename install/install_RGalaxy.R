#
# One off installation script to install RGalaxy on the server.
#

# Workaround for errors during RGalaxy install:
# Error : package ‘(name)’ was built before R 3.0.0: please re-install it
if (!require("stringr")) {
  install.packages("stringr", repos="http://cran.rstudio.com/") 
  library("stringr")
}
if (!require("digest")) {
  install.packages("digest", repos="http://cran.rstudio.com/") 
  library("digest")
}
if (!require("VennDiagram")) {
  install.packages("VennDiagram", repos="http://cran.rstudio.com/") 
  library("VennDiagram")
}

source("http://bioconductor.org/biocLite.R")
if(!require(RGalaxy)) {
  biocLite("RGalaxy")
}

if (!require("igraph")) {
  install.packages("igraph", repos="http://cran.rstudio.com/") 
  library("igraph")

}
if (!require("ggplot2")) {
  install.packages("ggplot2", repos="http://cran.rstudio.com/") 
  library("ggplot2")
}


# RGalaxy install problems:
# * Cannot find xml2-config ERROR: configuration failed for package ‘XML’
#    apt-get install libxml2-dev
