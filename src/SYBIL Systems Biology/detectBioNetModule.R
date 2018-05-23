##
## detect functional module by BioNet
##
suppressPackageStartupMessages(library(RGalaxy))

# TODO: Currently doesn't fit into the 'Downstream pathway analysis' workflow as that takes
# adjusted p-values but this function only works on non-adjusted pvalues.

detectBioNetModule <- function(
  differentialExpression = GalaxyInputFile(required=TRUE),
  interactome = GalaxyInputFile(required=TRUE),
  foldChangeOnly = GalaxyLogicalParam(),
  moduleNodes = GalaxyOutput("moduleNodes", "tabular")) {
  
  if (foldChangeOnly) {
    write.table('No p-values', file=moduleNodes, sep = '\t')
    stop()
  }
  
  differentialExpression <- read.table(differentialExpression, sep='\t', header = TRUE, as.is = TRUE)
  interactome <- read.table(interactome, sep='\t', as.is=TRUE)
  
  suppressPackageStartupMessages(library(igraph))
  suppressPackageStartupMessages(library(BioNet))
  
  g1=graph.data.frame(interactome ,directed=F)
  
  subnet<-subNetwork(differentialExpression[,1],g1)
  subnet<-rmSelfLoops(subnet)
    
  # To score each node of the network we fit a Beta-uniform mixture model
  # (BUM) [10] to the p-value distribution and subsequently use the parameters of
  # the model for the scoring function [5]. A false-discovery rate (FDR) of 0.00001 is chosen.
  pvals = differentialExpression[,3]
  names(pvals) = differentialExpression[,1]
  fb<-fitBumModel(pvals, plot=FALSE)
    
  # We set the False Discovery Rate(FDR)cutoff as 0.00001 other than the default value 0.001,in order to reduce the size of the subnetwork.  
  scores<-scoreNodes(subnet, fb, fdr=0.001)
  module<-runFastHeinz(subnet,scores)
  g<-getEdgeList(module)
  d<-as.matrix(g)
    
  BioNetModuleNode<-unique(c(d[,1],d[,2]))
  write.table(BioNetModuleNode, file=moduleNodes, row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
}