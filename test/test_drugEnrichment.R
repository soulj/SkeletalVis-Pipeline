source('../src/SYBIL Systems Biology/drugEnrichment.R')
sourcefile1 <- 'resources/RNASeq/mouseFoldChangeTable.txt'
sourcefile2<-"resources/Homology.mouse.txt"


drugs.mimic<- tempfile()
drugs.reverse<- tempfile()

drugEnrichment(differentialExpression = sourcefile1,species = "Mouse",homology = sourcefile2,foldChangeOnly = T,enrichedDrugsMimic = drugs.mimic,enrichedDrugsReverse = drugs.reverse)


drugs.mimic  <- read.delim(drugs.mimic)
head(drugs.mimic)

drugs.reverse  <- read.delim(drugs.reverse)
head(drugs.reverse)
