source('../src/SYBIL Systems Biology/GigaSubnetworks.R')
sourcefile1 <- 'resources/RNASeq/foldChangeTable.txt'
sourcefile2 <- 'resources/RNASeq/mouseStringNetwork.txt'
sourcefile3 <- 'resources/RNASeq/mouseGeneLengths.txt'

moduleNodes <- tempfile()
modulePlots <- tempfile()
visNetworks = tempfile()
summaryTable = tempfile()

GigaSubnetworks(differentialExpression = sourcefile1,interactome = sourcefile2,  
                species = "Mouse", foldChangeOnly = F, moduleNodes = moduleNodes, modulePlots = modulePlots,
                visNetworks = visNetworks, summaryTable = summaryTable)


sessionInfo()

moduleNodes  <- read.delim(moduleNodes)
dim(moduleNodes)
head(moduleNodes)

summaryTable  <- read.delim(summaryTable)
head(summaryTable)

visNetworks<-load(visNetworks)
get(visNetworks)[[1]]

