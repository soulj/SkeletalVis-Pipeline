suppressPackageStartupMessages(library(RGalaxy))

differentialExpressionPageRank <- function(
  differentialExpression =  GalaxyInputFile(required=TRUE),
  interactome =  GalaxyInputFile(required=TRUE),
  alpha = GalaxyNumericParam(0.7,required=F),
  backgroundSubtract =  GalaxyLogicalParam(checked=TRUE),
  rankedGenes  = GalaxyOutput("rankedGenes", "tabular") 
) {
  
  suppressPackageStartupMessages(library(igraph))
  suppressPackageStartupMessages(library(dplyr))
  
  #set up the networks
  interactome <- read.table(interactome, sep='\t', as.is=TRUE)
  interactome <- graph_from_data_frame(interactome ,directed=F)
  
  #filter the network based on expressed genes and extract the largest connected component
  differentialExpression<-read.delim(differentialExpression,stringsAsFactors = F)
  differentialExpression<-na.omit(differentialExpression)
  
  colnames(differentialExpression)[1]<-"GeneSymbol"
  differentialExpression<-  differentialExpression %>%
    group_by(GeneSymbol) %>%
    dplyr::slice(which.max(abs(2))) %>% as.data.frame
  
  presentList<-na.omit(match(differentialExpression[,1],V(interactome)$name))
  interactome<-induced.subgraph(interactome,presentList)
  interactome<-decompose.graph(interactome)
  interactome<-interactome[[which.max(sapply(interactome, vcount))]]
  
  #filter and order the expression table based on genes in the network
  presentList<-na.omit(match(V(interactome)$name,differentialExpression[,1]))
  differentialExpression<-differentialExpression[presentList,]
  
  #add the absolute fold change to the network
  V(interactome)$FC<-abs(as.numeric(differentialExpression[,2]))
  
  # calcuate the personalized PageRank score of all nodes in the network using the fold changes
  node.score <- stack(page.rank(interactome,directed = FALSE,damping=alpha, personalized =V(interactome)$FC,weights =E(interactome)$V3, options=NULL)$vector)

  if (backgroundSubtract == TRUE){
    # calcuate the background PageRank score of all nodes in the network
    node.score.background <- stack(page.rank(interactome,directed = FALSE,damping=alpha,weights =E(interactome)$V3, options=NULL)$vector)
    
    node.score$values <- node.score$values - node.score.background$values
    
  }

  node.score$rank<-rank(-node.score$values)
  colnames(node.score)<-c("score","gene","rank")
  
  #merge with the differential expression
  differentialExpression<-merge(differentialExpression,node.score,by.x=1,by.y="gene")
  n <- ncol(differentialExpression)
  differentialExpression<-differentialExpression[order(differentialExpression$rank,decreasing = F),]
  
  write.table(differentialExpression, file=rankedGenes , row.names=FALSE,quote=FALSE, sep="\t")
}