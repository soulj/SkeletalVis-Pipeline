source('../src/SYBIL Systems Biology/annotateSecretedProteins.R')
sourcefile1 <- 'resources/RNASeq/foldChangeTable.txt'
sourcefile2 <- 'resources/uniprot-secreted-mouse.txt'

secretedProteins <- tempfile()

annotateSecretedProteins(differentialExpression = sourcefile1,secretedReference = sourcefile2, secretedProteins =  secretedProteins)
