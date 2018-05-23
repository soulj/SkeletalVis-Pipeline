source('../src/SYBIL Systems Biology/differentialExpressionPageRank.R')
sourcefile1 <- 'resources/RNASeq/foldChangeTable.txt'
sourcefile2 <- 'resources/RNASeq/mouseStringNetwork.txt'
alpha <- 0.7

rankedGenes <- tempfile()

differentialExpressionPageRank(differentialExpression = sourcefile1,interactome = sourcefile2, alpha = alpha, backgroundSubtract = T,rankedGenes = rankedGenes)

head(read.delim(rankedGenes))