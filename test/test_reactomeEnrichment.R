source('../src/SYBIL Systems Biology/reactomeEnrichment.R')
sourcefile1 <- 'resources/RNASeq/mouseFoldChangeTable.txt'
sourcefile2 <- 'resources/RNASeq/mouseGeneLengths.txt'

pathways <- tempfile()
enrichmentPlot <- tempfile()

reactomeEnrichment(differentialExpression = sourcefile1,geneExonLengths = sourcefile2,species = "Mouse",foldChangeOnly = T,pathways = pathways,enrichmentPlot = enrichmentPlot)

sessionInfo()

pathways  <- read.delim(pathways)
print("pathways : ")
dim(pathways)
head(pathways)