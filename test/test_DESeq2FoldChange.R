source('../src/SYBIL Systems Biology/DESeq2FoldChange.R')
sourcefile1 <- 'resources/RNASeq/KallistoGeneAbundanceOutput.rdata'
sourcefile2 <- 'resources/RNASeq/sampleTable.txt'
sourcefile3 <- 'resources/RNASeq/comparisonTable.txt'

foldChangeTable <- tempfile()

DESeq2FoldChange(txiData = sourcefile1,
          sampleTable = sourcefile2,comparisonsTable = sourcefile3, foldChangeOnly = F,species = "Mouse",interaction = F,technicalReplicates=F,batchCorrect = F,supervised = F,foldChangeTable = foldChangeTable)

sessionInfo()

foldChangeTable  <- read.delim(foldChangeTable)
print("foldChangeTable : ")
dim(foldChangeTable)
head(foldChangeTable)