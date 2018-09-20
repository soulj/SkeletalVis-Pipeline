suppressPackageStartupMessages(library("RGalaxy"))

pathwaySlimEnrichment <- function(
  differentialExpression = GalaxyInputFile(required=TRUE),
  geneExonLengths = GalaxyInputFile(required = T, formatFilter = "tabular"),
  species=GalaxySelectParam(c("Human","Mouse","Rat","Horse","Zebrafish","Cow","Pig")),
  foldChangeOnly = GalaxyLogicalParam(),
  foldchange = GalaxyNumericParam(1.5),
  padj = GalaxyNumericParam(0.05),
  slimEnrichmentPathways     = GalaxyOutput("slimEnrichmentPathways", "tabular"),
  slimEnrichmentPlot     = GalaxyOutput("slimEnrichmentPlot", "png")) {
  
  suppressPackageStartupMessages(library("goseq"))
  suppressPackageStartupMessages(library("pathwaySlim"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("dplyr"))
  
  
  getEnrichedPathways<-function(data,sigData,geneExonLengths,species,threshold=0.05){
    
    data<-as.data.frame(data)
    sigData<-as.data.frame(sigData)
    
    genes<-ifelse(as.data.frame(data)[,1]%in% as.data.frame(sigData)[,1],1,0)
    names(genes)<-data[,1]
    
    #map from gene symbols to PATHWAY
    if (species == "Human") {
      pathway<-pathway2humangene
    }
    else if (species == "Mouse"){
      pathway<-pathway2mousegene
    }
    else if (species == "Zebrafish"){
      pathway<-pathway2zebrafishgene
    }
    else if (species == "Pig"){
      pathway<-pathway2piggene
    }
    else if (species == "Horse"){
      pathway<-pathway2horsegene
    } else if (species == "Cow"){
      pathway<-pathway2cowgene
    } else {
      pathway<-pathway2ratgene
    }
    
    geneExonLengths<-geneExonLengths[match(data[,1],geneExonLengths$gene_name),]
    
    x<-nullp(genes,bias.data = geneExonLengths$length,plot.fit=T)
    
    PATHWAY=goseq(x,gene2cat=pathway)
    PATHWAY$padj=p.adjust(PATHWAY$over_represented_pvalue,method="BH")
    PATHWAY$Enrichment=PATHWAY$numDEInCat/PATHWAY$numInCat*100
    PATHWAY$Adjpvaluelog=-log10(PATHWAY$padj)
    PATHWAY$Term<-unlist(PATHWAY$Term)
    PATHWAY.sig=PATHWAY[PATHWAY$padj<=threshold,]
    
    if (nrow(PATHWAY.sig)==0) {
      print("no significant pathways")
      return()
    }
    
    pathwayResults=list()
    
    for ( i in 1:nrow(PATHWAY.sig)) {
      
      #search pathway for the pathway term of interest and filter by differentially expressed genes
      pathwayTerm=PATHWAY.sig$category[i]
      termIDs=pathway[[pathwayTerm]]
      
      sig=sigData[sigData[,1] %in% termIDs ,]
      pathwayResults[[pathwayTerm]]=sig[,1]
    }
    names(pathwayResults)=PATHWAY.sig$category
    
    pathwayResults=lapply(pathwayResults,function(x) paste(x, sep="", collapse=" ") )
    pathwayResults=data.frame(Term= names(pathwayResults),Genes = unlist(pathwayResults),Adj.pvalue=PATHWAY.sig$padj,Enrichment=PATHWAY.sig$Enrichment,Adjpvaluelog=PATHWAY.sig$Adjpvaluelog)
    return(pathwayResults)
    
  }
  
  differentialExpression<-na.omit(read.delim(differentialExpression))
  geneExonLengths<-read.delim(geneExonLengths)
  
  colnames(differentialExpression)[1]<-"GeneSymbol"
  
  if (foldChangeOnly==TRUE){
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      slice(which.max(abs(2))) %>% as.data.frame
    
    
    differentialExpression.sig<-na.omit(differentialExpression[ abs(differentialExpression[,2])>=log2(foldchange),])
  } else{
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      slice(which.min(3)) %>% as.data.frame
    differentialExpression.sig<-na.omit(differentialExpression[ abs(differentialExpression[,2])>=log2(foldchange) & differentialExpression[,3] <=padj,])
  }
  
  #get the enriched pathways
  pathwayResults<-getEnrichedPathways(differentialExpression,differentialExpression.sig,geneExonLengths,species)
  
  #generate the enrichment plot
  g<-ggplot(pathwayResults,aes( x= Term, y= Enrichment,fill=Adjpvaluelog))+ geom_bar(stat ="identity")+ylab(label="Enrichment")+xlab(label="Pathway") + theme_classic(base_size=15) +  coord_flip() +  scale_x_discrete(limits = pathwayResults$Term)
  
  ggsave(filename = slimEnrichmentPlot, plot = g, device="png")
  
  write.table(pathwayResults,file=slimEnrichmentPathways,col.names=T,row.names=F,sep="\t",quote=F)
  
  
}
