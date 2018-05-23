source('../src/SYBIL Systems Biology/TF_Enrichment.R')
sourcefile1 <- 'resources/RNASeq/foldChangeTable.txt'

t<-tempfile()

TF_Enrichment(differentialExpression = sourcefile1,species = "Mouse",foldChangeOnly = T,enrichmentTable = t)

t<-read.delim(t)
dim(t)
head(t)
