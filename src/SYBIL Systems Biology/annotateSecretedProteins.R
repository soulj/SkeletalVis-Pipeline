suppressPackageStartupMessages(library("RGalaxy"))

annotateSecretedProteins <- function(
  differentialExpression = GalaxyInputFile(required=TRUE),
  secretedReference = GalaxyInputFile(required=TRUE),
  secretedProteins    = GalaxyOutput("secretedProteins", "tabular")
  
) {
  
  differentialExpression<-read.delim(differentialExpression)
  secretedReference <- read.delim( secretedReference,header=F)
  
  differentialExpression$secreted<-ifelse(differentialExpression[,1] %in%  secretedReference[,1],1,0)
  
  write.table(differentialExpression,file=secretedProteins,col.names=T,row.names=F,sep="\t",quote=F)
  
  
}
