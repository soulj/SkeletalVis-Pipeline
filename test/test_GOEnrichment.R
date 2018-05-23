source('../src/SYBIL Systems Biology/GOEnrichment.R')
sourcefile1 <- 'resources/RNASeq/mouseFoldChangeTable.txt'
sourcefile2 <- 'resources/RNASeq/mouseGeneLengths.txt'

GOTerms<- tempfile()
GOTerms.reduced<- tempfile()
mdsPlot<- tempfile()

GOEnrichment(differentialExpression = sourcefile1,geneLengths = sourcefile2,species = "Mouse",foldChangeOnly = T,enrichedTerms = GOTerms,enrichedTermsReduced = GOTerms.reduced,mdsPlot = mdsPlot)

sessionInfo()

GOTerms  <- read.delim(GOTerms)
head(GOTerms)

GOTerms.reduced  <- read.delim(GOTerms.reduced)
head(GOTerms.reduced)
