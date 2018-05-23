source('../src/SYBIL Systems Biology/PCARNASeq.R')
sourcefile1 <- 'resources/RNASeq/KallistoGeneAbundanceOutput.rdata'
sourcefile2 <- 'resources/RNASeq/sampleTable.txt'

pcaPlot <- tempfile()
geneInfluence <- tempfile()


PCARNASeq(txiData = sourcefile1,
          sampleTable = sourcefile2,species = "Mouse",pcaPlot = pcaPlot,technicalReplicates = F,batchCorrect = F,geneInfluence = geneInfluence)

geneInfluence <- read.delim(geneInfluence)
print("GeneInfluence: ")
dim(geneInfluence)
head(geneInfluence)
