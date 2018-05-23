##
## generate personalizedPageRank
##
suppressPackageStartupMessages(library(RGalaxy))

personalizedPageRank <- function(
  topicGenes   = GalaxyInputFile(required=TRUE),
  networkEdges = GalaxyInputFile(required=TRUE),
  rankedGenes  = GalaxyOutput("rankedGenes", "tabular") 
) {
  
  suppressPackageStartupMessages(library(igraph))
  
  interactome <- read.table(networkEdges, sep='\t', as.is=TRUE)
  g1=graph_from_data_frame(interactome ,directed=F)
  network.node<-V(g1)$name
  
  SignificantSeedGene <- read.table(topicGenes, sep='\t', as.is=TRUE)
  
  input<-as.character(SignificantSeedGene[,1])
  n<-length(network.node[network.node %in% input])
  
  
  ## set the value of personalized parameter, the sum of it is 1
  topicNode<-ifelse(network.node %in% input,1/n,0)
  ## calcuate the personalized PageRank score of all nodes in the network
  node.score<-page.rank (g1,algo="prpack",vids = V(g1), directed = FALSE,damping=0.85, personalized =topicNode,weights =NULL, options=NULL)
 
  dictionary<-as.matrix(node.score$vector)
  dictionary1<-cbind(rownames(dictionary),as.numeric(dictionary[,1]))
  
  #remove row with score of NaN, not sure if it is needed
  score1<-dictionary1[!is.nan(dictionary1[,1]),] 
  ## sort the score in a decreasing order
  sorted.score<-score1[order(as.numeric(score1[,2]),decreasing=TRUE),]
  write.table(sorted.score, file=rankedGenes , row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
}