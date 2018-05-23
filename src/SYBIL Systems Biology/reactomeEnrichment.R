suppressPackageStartupMessages(library("RGalaxy"))

reactomeEnrichment <- function(
  differentialExpression = GalaxyInputFile(required=TRUE),
  geneExonLengths = GalaxyInputFile(required = T, formatFilter = "tabular"),
  species=GalaxySelectParam(c("Human","Mouse")),
  foldChangeOnly = GalaxyLogicalParam(),
  padj = GalaxyNumericParam(0.05),
  pathways     = GalaxyOutput("pathways", "tabular"),
  enrichmentPlot     = GalaxyOutput("enrichmentPlot", "png")) {
  
  suppressPackageStartupMessages(library("goseq"))
  suppressPackageStartupMessages(library("reactome.db"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library("org.Mm.eg.db"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("dplyr"))
  

  
  getEnrichedReactomePathways<-function(data,sigData,geneExonLengths,species=c("Human","Mouse"),threshold=0.05){
    
    data<-as.data.frame(data)
    sigData<-as.data.frame(sigData)
    
    genes<-ifelse(as.data.frame(data)[,1]%in% as.data.frame(sigData)[,1],1,0)
    names(genes)<-data[,1]
    
    #map from gene symbols to REACTOME
    if (species == "Human") {
      eg2symbol=as.list(org.Hs.egSYMBOL2EG)
    }    else {
      eg2symbol=as.list(org.Mm.egSYMBOL2EG)
    }
    
    eg2reactome=as.list(reactomeEXTID2PATHID)
    grepREACTOME=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
    reactome=lapply(eg2symbol,grepREACTOME,eg2reactome)
    
    geneExonLengths<-geneExonLengths[match(data[,1],geneExonLengths$gene_name),]
    
    x<-nullp(genes,bias.data = geneExonLengths$length,plot.fit=T)
    
    REACTOME=goseq(x,gene2cat=reactome)
    REACTOME$padj=p.adjust(REACTOME$over_represented_pvalue,method="BH")
    xx <- as.list(reactomePATHID2NAME)
    REACTOME$Term=apply(REACTOME,1,function(x) xx[[unlist(x[1])]])
    REACTOME$Enrichment=REACTOME$numDEInCat/REACTOME$numInCat*100
    REACTOME$Adjpvaluelog=-log10(REACTOME$padj)
    REACTOME$Term[sapply(REACTOME$Term, length) == 0] <- NA
    REACTOME$Term<-unlist(REACTOME$Term)
    REACTOME.sig=REACTOME[REACTOME$padj<=threshold,]
    
    if (nrow(REACTOME.sig)==0) {
    print("no significant pathways")
    return()
    }
    
    reactomeResults=list()
    
    for ( i in 1:nrow(REACTOME.sig)) {
      
      #search reactome for the reactome term of interest and filter by differentially expressed genes
      reactomeTerm=REACTOME.sig$category[i]
      index=sapply(reactome,function(x) reactomeTerm %in% x)
      termIDs=names(index[index=="TRUE"])
      
      sig=sigData[sigData[,1] %in% termIDs ,]
      reactomeResults[[reactomeTerm]]=sig[,1]
    }
    names(reactomeResults)=REACTOME.sig$Term
    
    reactomeResults=lapply(reactomeResults,function(x) paste(x, sep="", collapse=" ") )
    reactomeResults=data.frame(Term= names(reactomeResults),Genes = unlist(reactomeResults),Adj.pvalue=REACTOME.sig$padj,Enrichment=REACTOME.sig$Enrichment,Adjpvaluelog=REACTOME.sig$Adjpvaluelog)
    reactomeResults[,1]<-gsub( "Homo sapiens: ", "", as.character(reactomeResults[,1]))
    
    return(reactomeResults)
    
  }
  
  differentialExpression<-read.delim(differentialExpression)
  geneExonLengths<-read.delim(geneExonLengths)
  colnames(differentialExpression)[1]<-"GeneSymbol"
  
  if (foldChangeOnly==TRUE){
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      slice(which.max(abs(2))) %>% as.data.frame
    
    
    differentialExpression.sig<-na.omit(differentialExpression[ abs(differentialExpression[,2])>=log2(1.5),])
  } else{
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      slice(which.min(3)) %>% as.data.frame
    differentialExpression.sig<-na.omit(differentialExpression[ abs(differentialExpression[,2])>=log2(1.5) & differentialExpression[,3] <=padj,])
  }
  
  #get the enriched pathways
  reactomeResults<-getEnrichedReactomePathways(differentialExpression,differentialExpression.sig,geneExonLengths,species)
  
  #generate the enrichment plot
  g<-ggplot(reactomeResults,aes( x= Term, y= Enrichment,fill=Adjpvaluelog))+ geom_bar(stat ="identity")+ylab(label="Enrichment")+xlab(label="Pathway") + theme_classic(base_size=15) +  coord_flip() +  scale_x_discrete(limits = reactomeResults$Term)
  
  ggsave(filename = enrichmentPlot, plot = g, device="png")
  
  write.table(reactomeResults,file=pathways,col.names=T,row.names=F,sep="\t",quote=F)
  
  
}
