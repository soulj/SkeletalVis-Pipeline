#
# One off installation script for BioConductor libraries. Several tools will require this.
# Run as root.
# Java needs to be correctly configred using:
# R CMD javareconf
#

# ensure the rJava package is available too
# install.packages("rJava")
# and igraph
# install.packages("igraph")

# Upgrading BioConductor:
#  One way of achieving this is to remove the current version,
#  allowing biocLite() to then install the latest.
# remove.packages("BiocInstaller")

source("http://bioconductor.org/biocLite.R")

# On Debian, may also need to install these packages to compile:
# libblas-dev
# liblapack-dev

# If rsqlite fails to compile then you can also install package:
# r-cran-rsqlite

# Install basic bioconductor packages and update any already installed.
biocLite(ask = FALSE)

if(!require(affy)) {
  biocLite("affy")
  # CDF is required to complete the background-correction step.
  # Note that affy will attempt to detect and load the correct CDF if available, but
  # installation will fail unless you are root. An alternative is to install the CDF required manually.
  biocLite("mouse4302cdf")
}

# AnnotationDbi also required to complete the background-correction step.
if(!require(AnnotationDbi)) {
  biocLite("AnnotationDbi")
}

# For identifyDifferentialGeneExpressions
if(!require(lumi)) {
  biocLite("lumi")
}
if(!require(lumiHumanIDMapping)) {
  biocLite("lumiHumanIDMapping")
}
if(!require(lumiHumanAll.db)) {
  biocLite("lumiHumanAll.db")
}
if(!require(limma)) {
  biocLite("limma")
}
if(!require(gcrma)) {
  biocLite("gcrma")
}
if(!require(affyPLM)) {
  biocLite("affyPLM")
}
if(!require(affyQCReport)) {
  biocLite("affyQCReport")
}
if(!require(GenomicFeatures)) {
  biocLite("GenomicFeatures")
}
if(!require(hgu95av2.db)) {
  biocLite("hgu95av2.db")
}
if(!require(mouse430a2.db)) {
  biocLite("mouse430a2.db")
}
if(!require(mouse4302.db)) {
  biocLite("mouse4302.db")
}
if(!require(mouse430a2cdf)) {
  biocLite("mouse430a2cdf")
}
if(!require(moe430acdf)) {
  biocLite("moe430acdf")
}
if(!require(STRINGdb)) {
  biocLite("STRINGdb")
}
if(!require(GO.db)) {
  biocLite("GO.db")
}
if(!require(paxtoolsr)) {
  biocLite("paxtoolsr")
}
if(!require(biomaRt)) {
  biocLite("biomaRt")
}
if(!require(igraph)) {
  biocLite("igraph")
}
if(!require(graph)) {
  biocLite("graph")
}
if(!require(RBGL)) {
  biocLite("RBGL")
}
if(!require(BioNet)) {
  biocLite("BioNet")
}
if(!require(GO.db)) {
  biocLite("GO.db")
}
if(!require(topGO)) {
  biocLite("topGO")
}
if(!require(clusterProfiler)) {
  biocLite("clusterProfiler")
}
if(!require(DESeq2)) {
  biocLite("DESeq2")
}
if(!require(EnsDb.Hsapiens.v79)) {
  biocLite("EnsDb.Hsapiens.v79")
}
if(!require(EnsDb.Mmusculus.v79)) {
  biocLite("EnsDb.Mmusculus.v79")
}
if(!require(goseq)) {
  biocLite("goseq")
}
if(!require(genefilter)) {
  biocLite("genefilter")
}
if(!require(org.Hs.eg.db)) {
  biocLite("org.Hs.eg.db")
}
if(!require(org.Mm.eg.db)) {
  biocLite("org.Mm.eg.db")
}
if(!require(tximport)) {
  biocLite("tximport")
}
if(!require(reactome.db)) {
  biocLite("reactome.db")
}
if(!require(GEOquery)) {
  biocLite("GEOquery")
}
if(!require(oligo)) {
  biocLite("oligo")
}
if(!require(vsn)) {
  biocLite("vsn")
}
if(!require(lumiMouseIDMapping)) {
  biocLite("lumiMouseIDMapping")
}
if(!require(lumiHumanIDMapping)) {
  biocLite("lumiHumanIDMapping")
}
if(!require(mgug4122a.db)) {
  biocLite("mgug4122a.db")
}
if(!require(hgu133plus2.db)) {
  biocLite("hgu133plus2.db")
}
#install the libraries needed to install RcisTarget - a developmental package
if(!require("devtools")){
  install.packages("devtools", repos="http://cran.rstudio.com/")
}
if(!require("zoo")){
  install.packages("zoo", repos="http://cran.rstudio.com/")
}
if(!require("AUCell")){
  biocLite("AUCell")
}
if(!require("RcisTarget")){
  devtools::install_github("aertslab/RcisTarget")
}
#install the required annotation packages directly from the lab's website
if(!require("RcisTarget.mm9.motifDatabases.20k")) {
  install.packages("http://scenic.aertslab.org/downloads/databases/RcisTarget.mm9.motifDatabases.20k_0.1.1.tar.gz",repos = NULL,type="source")
}
if(!require("RcisTarget.hg19.motifDatabases.20k")) {
  install.packages("http://scenic.aertslab.org/downloads/databases/RcisTarget.hg19.motifDatabases.20k_0.1.1.tar.gz",repos = NULL,type="source")
}
if(!require("cowplot")){
  install.packages("cowplot", repos="http://cran.rstudio.com/")
}
if(!require("squash")){
  install.packages("squash", repos="http://cran.rstudio.com/")
}
if(!require("intergraph")){
  install.packages("intergraph", repos="http://cran.rstudio.com/")
}
if(!require("ggnet")){
  devtools::install_github("briatte/ggnet")
}
if(!require("ggpubr")){
  install.packages("ggpubr", repos="http://cran.rstudio.com/")
}
if(!require("visNetwork")){
  install.packages("visNetwork", repos="http://cran.rstudio.com/")
}
if(!require("plotly")){
  install.packages("plotly", repos="http://cran.rstudio.com/")
}
if(!require("dplyr")){
  install.packages("dplyr", repos="http://cran.rstudio.com/")
}
if(!require("tidyr")){
  install.packages("tidyr", repos="http://cran.rstudio.com/")
}
if(!require("httr")){
  install.packages("httr", repos="http://cran.rstudio.com/")
}
if(!require("GO.db")){
  biocLite("GO.db")
}
if(!require("GOSemSim")){
  biocLite("GOSemSim")
}
if(!require("readr")){
  install.packages("readr", repos="http://cran.rstudio.com/")
}
if(!require("htmlwidgets")){
  install.packages("htmlwidgets", repos="http://cran.rstudio.com/")
}
if(!require("sva")){
  biocLite("sva")
}
if(!require("org.Rn.eg.db")){
  biocLite("org.Rn.eg.db")
}
if(!require("mogene10sttranscriptcluster.db")){
  biocLite("mogene10sttranscriptcluster.db")
}
if(!require("arrayQualityMetrics")){
  biocLite("arrayQualityMetrics")
}
if(!require("ENAbrowseR")){
  devtools::install_github("cstubben/ENAbrowseR")
}

if(!require("rhdf5")){
  biocLite("rhdf5")
}
if(!require("GeoDE")){
  biocLite("GeoDE")
}
if(!require("org.Dr.eg.db")){
biocLite("org.Dr.eg.db")
}
if(!require("org.Rn.eg.db")){
  biocLite("org.Rn.eg.db")
}
if(!require("org.Bt.eg.db")){
  biocLite("org.Bt.eg.db")
}
if(!require("org.Ss.eg.db")){
  biocLite("org.Ss.eg.db")
}

#c("pd.hugene.1.0.st.v1","pd.hugene.2.0.st","pd.mogene.1.0.st.v1","pd.mogene.2.0.st","pd.ragene.1.0.st.v1","pd.ragene.2.0.st")

