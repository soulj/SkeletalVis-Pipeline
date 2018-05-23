options("width"=200)
source('../src/SYBIL Systems Biology/celltypeSpecificNetwork.R')
sourcefile1 <- 'resources/sampleGeneList.txt'
sourcefile2 <- 'resources/SampleInteractome.txt'
## not yet implemented
##quit()

edges <- tempfile()
isolatednodes <- tempfile()
celltypeSpecificNetwork(
  expressedGenes       = sourcefile1,
  interactome          = sourcefile2,
  networkEdges         = edges,
  networkIsolatedNodes = isolatednodes)

edgeTable <- read.table(edges, header=FALSE, fill=TRUE)
print("Edges:")
dim(edgeTable)
head(edgeTable)
isolatedNodesTable <- read.table(isolatednodes, header=FALSE, fill=TRUE)
print("Isolated nodes: ")
dim(isolatedNodesTable)
head(isolatedNodesTable)
