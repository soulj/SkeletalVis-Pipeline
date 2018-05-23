##
## Create Consensus DB interactome
##
## This function is to retrieve protein-protein interactions from ConsensusDB.
## The input data is a protein-protein interaction file downloaded from ConsensusPathDB.
## Interactions in the input data may be of 3 different types as shown below:

##  BIND      14988405        5HT2C_MOUSE,DLG3_MOUSE  0.0149385
##  BIND      15603741        ALEX_MOUSE.GNAS1_MOUSE.GNAS2_MOUSE.GNAS3_MOUSE,ZDHC3_MOUSE      NA
##  BIND      12614612        ORC5_MOUSE      0.0089412599999999995

## When the participant node (ALEX_MOUSE.GNAS1_MOUSE.GNAS2_MOUSE.GNAS3_MOUSE) is a complex formed by several proteins, 
## then the confidence score is "NA" and the interaction is deleted.
## Interactions with only one node to represent selfregulation are extended to an interaction between two identical nodes.
## The names of participant nodes should be formated in the same way as the gene names from STRING by removing the end component, for example "_MOUSE".
## The expected output of the above three recordes should be:
##  5HT2C  DLG3
##  ORC5   ORC5
suppressPackageStartupMessages(library(RGalaxy))

createConsensusDBInteractome <- function(
  consensusDBData = GalaxyInputFile(),
  interactome     = GalaxyOutput("interactome", "tabular")
) {
  data <- read.table(consensusDBData, sep='\t', as.is=TRUE)
  # filter out rows with no confidence score 
  data <- data[!is.na(as.numeric(data$V4)),]
  # only interested in column 3
  data <- data[c(3)]
  # split on comma
  data<-read.csv(textConnection(data[["V3"]]), header=FALSE)
  # remove the species from the gene name, pattern is GENE_SPECIES
  data$V1<-gsub('_.*$', '', data$V1)
  data$V2<-gsub('_.*$', '', data$V2)
  #if the second column is empty, copy the first column into the second column
  data$V2<-ifelse(data$V2 == '', data$V1, data$V2)
  write.table(data, file=interactome, row.names=FALSE, col.names=FALSE, quote=FALSE, na="", sep="\t", append=TRUE)
}
