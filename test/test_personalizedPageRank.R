options("width"=200)
source('../src/SYBIL Systems Biology/personalizedPageRank.R')
sourcefile1 <- 'resources/PageRankTestTopicGenes.txt'
sourcefile2 <- 'resources/PageRankTestMLIIedges.txt'
## not yet implemented
##quit()

edges <- tempfile()

personalizedPageRank(
  topicGenes       = sourcefile1,
  networkEdges          = sourcefile2,
  rankedGenes          = edges
  )

rankedGenes  <- read.table(edges, header=FALSE, fill=TRUE)
print("Edges:")
dim(rankedGenes)
head(rankedGenes)

