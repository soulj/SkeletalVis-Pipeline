suppressPackageStartupMessages(library("RGalaxy"))

GOEnrichment <- function(
  differentialExpression = GalaxyInputFile(required=TRUE),
  geneLengths            = GalaxyInputFile(required=TRUE),
  foldChangeOnly         = GalaxyLogicalParam(),
  species=GalaxySelectParam(c("Human","Mouse","Rat","Horse","Zebrafish","Cow","Pig")),
  foldchange = GalaxyNumericParam(1.5),
  padj = GalaxyNumericParam(0.05),
  enrichedTerms    = GalaxyOutput("enrichedTerms", "tabular"),
  enrichedTermsReduced   = GalaxyOutput("enrichedTerms.reduced", "tabular"),
  mdsPlot   = GalaxyOutput("GO.MDS", "html")
  
) {
  
  suppressPackageStartupMessages(library("goseq"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library("org.Mm.eg.db"))
  suppressPackageStartupMessages(library("org.Rn.eg.db"))
  suppressPackageStartupMessages(library("org.Dr.eg.db"))
  suppressPackageStartupMessages(library("org.Ss.eg.db"))
  suppressPackageStartupMessages(library("org.Bt.eg.db"))
  suppressPackageStartupMessages(library("org.Ecaballus.eg.db"))
  suppressPackageStartupMessages(library("GO.db"))
  suppressPackageStartupMessages(library("GOSemSim"))
  suppressPackageStartupMessages(library("tidyr"))
  suppressPackageStartupMessages(library("plotly"))
  suppressPackageStartupMessages(library("htmlwidgets"))
  suppressPackageStartupMessages(library("dplyr"))
  
  getEnrichedGOTerms<-function(diffExp,sigDiffExp,geneLengths,species,threshold=0.05){
    
    #map from gene symbols to REACTOME
    if (species == "Human") {
      gene2GO<-AnnotationDbi::select(org.Hs.eg.db,keys(org.Hs.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Mouse"){
      gene2GO<-AnnotationDbi::select(org.Mm.eg.db,keys(org.Mm.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Cow"){
      gene2GO<-AnnotationDbi::select(org.Bt.eg.db,keys(org.Bt.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Horse"){
      gene2GO<-AnnotationDbi::select(org.Ecaballus.eg.db,keys(org.Ecaballus.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Zebrafish"){
      gene2GO<-AnnotationDbi::select(org.Dr.eg.db,keys(org.Dr.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Pig"){
      gene2GO<-AnnotationDbi::select(org.Ss.eg.db,keys(org.Ss.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else{
      gene2GO<-AnnotationDbi::select(org.Rn.eg.db,keys(org.Rn.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    }
    
    gene2GO<-unstack(gene2GO[,c(1,5)])
    
    genes<-ifelse(diffExp[,1] %in% sigDiffExp[,1] ,1,0)
    names(genes)<-diffExp[,1]
    
    geneLengths<-geneLengths[match(differentialExpression[,1],geneLengths$gene_name),]
    x<-nullp(genes,bias.data = geneLengths$length,plot.fit=F)
    GOTerms <- goseq(x,gene2cat=gene2GO)
    
    #filter to keep just BP GOTerms
    GOTerms <- GOTerms[ GOTerms$ontology =="BP",]
    GOTerms$padj<-p.adjust(GOTerms$over_represented_pvalue,method="BH")
    
    GOTerms.sig <- GOTerms[GOTerms$padj<=threshold,]
    
    
    if (nrow(GOTerms.sig)==0) {
      print("no significant GO terms!")
      return()
    }
    GOTerms.sig$enrichment <- GOTerms.sig$numDEInCat/GOTerms.sig$numInCat
    
    GOResults=list()
    
    for ( i in 1:nrow(GOTerms.sig)) {
      
      #search reactome for the reactome term of interest and filter by differentially expressed genes
      GOTerm <- GOTerms.sig$category[i]
      index <- sapply(gene2GO ,function(x) GOTerm %in% x)
      termIDs <- names(index[index=="TRUE"])
      
      sig <- sigDiffExp[sigDiffExp[,1] %in% termIDs ,]
      GOResults[[GOTerm]]=sig[,1]
    }
    names(GOResults)=GOTerms.sig$term
    
    GOResults <- lapply(GOResults,function(x) paste(x, sep="", collapse=" ") )
    GOResults <- data.frame(Term= names(GOResults),ID=GOTerms.sig$category,Genes = unlist(GOResults),Adj.pvalue=GOTerms.sig$padj,Enrichment=GOTerms.sig$enrichment)
    
    GOResults.reduced <-try(simplify(GOResults,gene2GO))
    
    return(list(GOResults=GOResults,GOResults.reduced=GOResults.reduced))
  }
  
  simplify <- function(GORes,gene2GO){
    
    if (species == "Human") {
      semData<-godata("org.Hs.eg.db","SYMBOL","BP")
    } else if (species == "Mouse") {
      semData<-godata("org.Mm.eg.db","SYMBOL","BP")
    } else if (species == "Zebrafish") {
      semData<-godata("org.Dr.eg.db","SYMBOL","BP")
    } else if (species == "Horse") {
      semData<-godata("org.Mm.eg.db","SYMBOL","BP")
    } else if (species == "Cow") {
      semData<-godata("org.Bt.eg.db","SYMBOL","BP")
    } else if (species == "Pig") {
      semData<-godata("org.Ss.eg.db","SYMBOL","BP")
    } else {
      semData<-godata("org.Rn.eg.db","SYMBOL","BP")
    }
    
    sim <- mgoSim(GORes$ID, GORes$ID,
                  semData = semData,
                  measure="Rel",
                  combine=NULL)
    
    sim[ is.na(sim)] <- 0
    
    ## to satisfy codetools for calling gather
    go1 <- go2 <- similarity <- NULL
    
    sim.df <- as.data.frame(sim)
    sim.df$go1 <- row.names(sim.df)
    sim.df <- gather(sim.df, go2, similarity, -go1)
    sim.df <- sim.df[!is.na(sim.df$similarity),]
    sim.df <- sim.df[ order(sim.df$similarity,decreasing = T),]
    
    #get the simiar term pairs
    sim.df <- sim.df[sim.df$go1 !=sim.df$go2,]
    sim.df <- sim.df[sim.df$similarity > 0.4,]
    
    #mark high fequency terms
    GO2Gene<-unstack(stack(gene2GO)[2:1])
    freq<-sapply(GO2Gene,length)
    freqCutOff<-length(gene2GO)*0.05
    highFreqTerms<-names(freq[freq>freqCutOff])
    
    sim.df$remove <- apply(sim.df,1,function(x) {
      if (x[1] %in% highFreqTerms){
        return(x[1])
      }
      if (x[2] %in% highFreqTerms){
        return(x[2])
      } else {
        return(NA)
        }
      })
    
    remove<-na.omit(sim.df$remove)
    sim.df<-sim.df[ is.na(sim.df$remove),]
    
    #iterate and remove the term with the least significant p-value
    sim.df$go1.pval<-GORes$Adj.pvalue[match(sim.df$go1,GORes$ID)]
    sim.df$go2.pval<-GORes$Adj.pvalue[match(sim.df$go2,GORes$ID)]
    
    childTerms<-as.list(GOBPCHILDREN)
    
    for (i in 1:nrow(sim.df)){
      
      #check to see if the goterm has already been marked for removal
      if (sim.df[i,"go1"] %in% remove){
        next
      }
      if (sim.df[i,"go2"] %in% remove){
        next
      }
      
      go1.pval <- sim.df[i,"go1.pval"]
      go2.pval <- sim.df[i,"go2.pval"]
      
      #if the p-values are equal then check if parent-child and keep the child term if so
      if (go1.pval==go2.pval){
        go1 <- sim.df[i,"go1"]
        go2 <- sim.df[i,"go2"]
        if(go2 %in% childTerms[[go1]]){
          remove<-c(remove,go2)
          next
        } else if (go1 %in% childTerms[[go2]])
          remove<-c(remove,go1)
          next
        }
      
      #remove least sig term
      remove<-c(remove,sim.df[i,which.max( c(go1.pval,go2.pval))])
      
      
    }

    GORes.filt <- GORes[ !GORes$ID %in% remove,]
    sim.filt <- sim[as.character(GORes.filt$ID),as.character(GORes.filt$ID)]
    
    fit <- cmdscale(1-sim.filt,eig=TRUE, k=2)
    x <- fit$points[,1]
    y <- fit$points[,2]
    
    GORes.filt.plot<-GORes.filt
    GORes.filt.plot$x <- x
    GORes.filt.plot$y <- y
    GORes.filt.plot$log10Adjpvalue<- -log10(GORes.filt.plot$Adj.pvalue)
    
    GO.MDS<-plot_ly(GORes.filt.plot, x = GORes.filt.plot$x, y = GORes.filt.plot$y,
            mode = "markers", type = "scatter",color = GORes.filt.plot$log10Adjpvalue,size=log2(GORes.filt.plot$Enrichment), text = paste("Term: ", GORes.filt$Term,""),marker = list(sizeref = 0.05))  %>%  colorbar( title='-log10 Adj.pvalue')
 
    return(list(GORes.filt,GO.MDS))
    
    }
  
  
  
    differentialExpression<-na.omit(read.delim(differentialExpression))
    geneLengths<-read.delim(geneLengths)

    
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
  GOResultsList<- getEnrichedGOTerms(differentialExpression,differentialExpression.sig,geneLengths,species)
  GOResults<-GOResultsList[[1]]
  write.table(GOResults,file=enrichedTerms,col.names=T,row.names=F,sep="\t",quote=F)
  
  if (!inherits(GOResultsList[[2]], "try-error")){
    write.table(GOResultsList[[2]][[1]],file=enrichedTermsReduced,col.names=T,row.names=F,sep="\t",quote=F)
    GO.MDS <- GOResultsList[[2]][[2]]
    GO.MDS <- saveWidget(GO.MDS,file=mdsPlot)
  }

}
