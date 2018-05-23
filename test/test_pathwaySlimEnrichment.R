source('../src/SYBIL Systems Biology/pathwaySlimEnrichment.R')
sourcefile1 <- 'resources/RNASeq/mouseFoldChangeTable.txt'
sourcefile2 <- 'resources/RNASeq/mouseGeneLengths.txt'

pathways <- tempfile()
enrichmentPlot <- tempfile()

pathwaySlimEnrichment (differentialExpression = sourcefile1,geneExonLengths = sourcefile2,species = "Mouse",foldChangeOnly = T,slimEnrichmentPathways = pathways,slimEnrichmentPlot = enrichmentPlot)

pathways  <- read.delim(pathways)
print("pathways : ")
dim(pathways)
head(pathways)