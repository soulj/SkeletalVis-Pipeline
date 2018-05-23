##
## generate celltypeSpecificNetwork
##
suppressPackageStartupMessages(library(RGalaxy))

celltypeSpecificNetwork <- function(
  expressedGenes       = GalaxyInputFile(required=TRUE),
  interactome          = GalaxyInputFile(required=TRUE),
  networkEdges         = GalaxyOutput("networkEdges", "tabular"),
  networkIsolatedNodes = GalaxyOutput("networkIsolatedNodes", "tabular")
) {
  
  ## This function will generate cell type specific network based on an input gene list and a two column interactome.
  ## The format of input gene list contained one gene symbol per row as below:
  
  ##  5HT2C
  ##  MDM2
  ##  DHC24
  ##  ORC5
  ##  CDT1
  ##  GEMI
  ##  RNF6
  
  ## The format of interactome contained two columns of gene symbol as below:
  
  ##  5HT2C	DLG3
  ##  MDM2	DHC24
  ##  ORC5	ORC5
  ##  CDT1	GEMI
  ##  ORC3	ORC3
  ##  HVM53	KV2A7
  ##  SHIP1	Q9D031
  ##  RNF6	LIMK1
  
  ## The function will produce two output: the edge list of network and a list isolated nodes in the network.
  ## The edge list is a two column interactome in which all nodes are found in the input gene list and connected by edges.
  ## The isolated nodes contained nodes found in the input gene list but they only interact with nodes absent in the input gene list.
  
  ## read the input interactome
  interactome <- read.table(interactome, sep='\t', as.is=TRUE)
  ## read the input gene list
  inputGene <- read.table(expressedGenes, sep='\t', as.is=TRUE)
  ## filter duplicate edges in the interactome
  interactome1= interactome[!duplicated(interactome),]
  ## get all nodes in the interactome
  allNode<-unique(c(interactome1[,1],interactome1[,2]))
  ## get all nodes in the cell type specific network
  networkNode<-intersect(inputGene[,1],allNode)
  
  ## get cell type specific network edges
  interactome2<-subset(interactome1,(interactome1[,1]%in%inputGene[,1])&(interactome1[,2]%in%inputGene[,1]))
  ## The output should be two columns of gene symbol.
  ## The networkEdges output of above example is:
  ##  MDM2	DHC24
  ##  ORC5	ORC5
  ##  CDT1	GEMI
  write.table(interactome2, file=networkEdges, row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
  ## get nodes present in the edges of cell type specific network
  edgeNode<-unique(c(interactome2[,1],interactome2[,2]))
  ## get network isolated nodes
  node2<-setdiff(networkNode,edgeNode)
  ## The  output should be one column of gene symbol
  ## The networkIsolatedNodes output of above example is:
  ##  5HT2C
  ##  RNF6
  write.table(node2, file=networkIsolatedNodes, row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
  }