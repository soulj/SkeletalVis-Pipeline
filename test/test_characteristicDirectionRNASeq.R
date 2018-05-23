source('../src/SYBIL Systems Biology/characteristicDirectionRNASeq.R')
sourcefile1 <- 'resources/RNASeq/KallistoGeneAbundanceOutput.rdata'
sourcefile2 <- 'resources/RNASeq/sampleTable.txt'
sourcefile3 <- 'resources/RNASeq/comparisonTable.txt'

foldChangeTable <- tempfile()

characteristicDirectionRNASeq(txiData = sourcefile1,
                 sampleTable = sourcefile2,comparisonsTable = sourcefile3, foldChangeOnly = F,species = "Mouse",technicalReplicates=F,batchCorrect = F,supervised = F,chrDirTable  = foldChangeTable)

sessionInfo()

foldChangeTable  <- read.delim(foldChangeTable)
print("foldChangeTable : ")
dim(foldChangeTable)
head(foldChangeTable)
