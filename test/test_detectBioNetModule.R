source('../src/SYBIL Systems Biology/detectBioNetModule.R')
sourcefile1 <- 'resources/differentialExpressionPValues.tabular'
sourcefile2 <- 'resources/RNASeq/MouseNetwork.txt'
modulenodes <- tempfile()

detectBioNetModule(
  differentialExpression = sourcefile1,
  interactome = sourcefile2,
  foldChangeOnly = FALSE,
  moduleNodes = GalaxyOutput("moduleNodes", "tabular"))

modulenodesTable <- read.table(modulenodes, header=FALSE, fill=TRUE)
dim(modulenodesTable)
head(modulenodesTable)

detectBioNetModule(
  differentialExpression = sourcefile1,
  interactome = sourcefile2,
  foldChangeOnly = TRUE,
  moduleNodes = GalaxyOutput("moduleNodes", "tabular"))
  
modulenodesTable <- read.table(modulenodes, header=FALSE, fill=TRUE)
dim(modulenodesTable)
head(modulenodesTable)