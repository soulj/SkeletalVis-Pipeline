suppressPackageStartupMessages(library("RGalaxy"))


DESeq2FoldChange<-function(txiData=GalaxyInputFile(required=T,formatFilter="rdata"),sampleTable=GalaxyInputFile(required=T,formatFilter="tabular"),comparisonsTable=GalaxyInputFile(required=T,formatFilter="tabular"),species=GalaxySelectParam(c("Human","Mouse","Rat","Pig","Cow","Horse","Zebrafish"),required=T),foldChangeOnly=GalaxyLogicalParam(),interaction=GalaxyLogicalParam(),technicalReplicates=GalaxyLogicalParam(),batchCorrect=GalaxyLogicalParam(),supervised=GalaxyLogicalParam(),foldChangeTable=GalaxyOutput("foldChangeTable","tabular")){
  
  
  suppressPackageStartupMessages(library("DESeq2"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("genefilter"))
  suppressPackageStartupMessages(library("EnsDb.Hsapiens.v79"))
  suppressPackageStartupMessages(library("EnsDb.Mmusculus.v79"))
  suppressPackageStartupMessages(library("EnsDb.Rnorvegicus.v79"))
  suppressPackageStartupMessages(library("EnsDb.Ecaballus.v79"))
  suppressPackageStartupMessages(library("EnsDb.Drerio.v79"))
  suppressPackageStartupMessages(library("EnsDb.Btaurus.v79"))
  suppressPackageStartupMessages(library("EnsDb.Sscrofa.v79"))
  
  suppressPackageStartupMessages(library("sva"))
  
  
  getFoldChangeDataFrame<-function(rld,sampleTable,condition,numerator,denominator) {
    
    
    print(condition)
    #get the samples which match to the numers and denoms
    numerator.ind<-which(sampleTable[,condition] %in% numerator)
    denominator.ind<-which(sampleTable[,condition] %in% denominator)
    
    print(numerator.ind)
    print(denominator.ind)
    
    rld<-assay(rld)
    
    #get the mean expression values for each gene if neccessary
    if (length(numerator.ind)>1){
      numerator.dat<-rowMeans(rld[,numerator.ind])
    } else {
      numerator.dat<-rld[,numerator.ind]
    }
    if (length(denominator.ind)>1){
      denominator.dat<-rowMeans(rld[,denominator.ind])
    } else {
      denominator.dat<-rld[,denominator.ind]
    }
    
    #calculate the fold change
    data<-data.frame(log2FC = numerator.dat -denominator.dat)
    colnames(data)<-paste(numerator,denominator,sep="vs")

    return(data)
    
  }
  
  
  getResultsDataFrame<-function(dds,condition,numerator,denominator) {
    
    data<-as.data.frame(lfcShrink(dds,contrast = c(condition,numerator,denominator)))
    data<-data[,c("log2FoldChange","padj")]
    colnames(data)<-paste(paste(numerator,denominator,sep="vs"),colnames(data),sep="_")
    return(data)
    
  }
  

  
  #load the input files
  load(txiData)
  sampleTable<-read.delim(sampleTable)
  
  #remove the file info
  sampleTable<-sampleTable[,colnames(sampleTable)!="File"]
  comparisonsTable<-read.delim(comparisonsTable,stringsAsFactors = F)
  
  #match the txi dataframes to the order of the sampleTable sample names
  txi$counts <- txi$counts[ ,match(as.character(sampleTable[,1]),colnames(txi$counts))]
  txi$abundance <- txi$abundance[ ,match(as.character(sampleTable[,1]),colnames(txi$abundance))]
  txi$length <- txi$length[ ,match(as.character(sampleTable[,1]),colnames(txi$length))]
  
  #get the design formula from the columns in the sampleTable (using the same order)
  factors<-colnames(sampleTable)[-1]
  

  
  

  #add an interaction term to the design formula if selected
  if (interaction){
    mainFactors<-paste(factors, collapse=" + ")
    interactionTerm<-paste(factors, collapse=":")
    designFormula <- as.formula(paste("~", paste(mainFactors, interactionTerm, sep=" + ")))
    for (i in 1:nrow(comparisonsTable)){
      factor<-as.character(comparisonsTable[i,"Factor"])
      baseLevel<-as.character(comparisonsTable[i,"baseLevel"])
      sampleTable[,factor]<-relevel(sampleTable[,factor],baseLevel)
    }
  } else {
    designFormula <- as.formula(paste("~", paste(factors, collapse=" + ")))
    
  }
  
  #perform normalisation
  
  #if technical replicates need collapsing
  if (technicalReplicates){
    dds<-DESeqDataSetFromTximport(txi = txi,colData = sampleTable,design = designFormula,ignoreRank = T)
    dds<-collapseReplicates(dds,paste0(sampleTable[,factors[1]],sampleTable[,factors[2]]))
    design(dds)<-as.formula(paste("~", paste(factors[-1], collapse=" + ")))
  } else {
    dds<-DESeqDataSetFromTximport(txi = txi,colData = sampleTable,design = designFormula)
  }
  
  #if only fold changes needed
  if (foldChangeOnly == TRUE) {
    
    rld <- rlogTransformation( dds )
    
    #make the fold change comparisons
    resultsTable<- lapply(1:nrow(comparisonsTable),function(x) getFoldChangeDataFrame(rld,sampleTable,comparisonsTable[x,1],comparisonsTable[x,2],comparisonsTable[x,3]))
    
  } else {
    
    
    
    if (interaction){
      dds<-DESeq(dds)
      resultsTable<-getInteractionResultsDataFrameList(dds)
      
    } else{
      
      if (batchCorrect){
        idx <- rowMeans(counts(dds)) > 1
        dat <- counts(dds)[idx,]
        mod <- model.matrix(designFormula, colData(dds))
        mod0 <- model.matrix(~ 1, colData(dds))
        
        if (supervised){
          dds<-DESeq(dds[idx,])
          resultsTable<-as.data.frame(results(dds,contrast = c(comparisonsTable[1,1],comparisonsTable[1,2],comparisonsTable[1,3])))
          resultsTable <- resultsTable[order(abs(resultsTable$pvalue)),]
          #take all but the top 5000 genes ranked by differential expression
          controls <-!rownames(dds) %in% rownames (resultsTable)[1:5000]
          
          svseq <- svaseq(dat, mod, mod0,controls=controls)
          
        } else{
        
        svseq <- svaseq(dat, mod, mod0)
        }
        
        colnames(svseq$sv)<-paste0("batch",1:svseq$n.sv)
        
        sampleTable <- cbind(svseq$sv,sampleTable)
        designFormula <- as.formula(paste("~", colnames(svseq$sv),"+" , paste(factors, collapse=" + ")))
        dds<-DESeqDataSetFromTximport(txi = txi,colData = sampleTable,design = designFormula)
      }
      dds<-DESeq(dds)
    
      #make the statistical comparisons
      resultsTable<- lapply(1:nrow(comparisonsTable),function(x) getResultsDataFrame(dds,comparisonsTable[x,1],comparisonsTable[x,2],comparisonsTable[x,3]))
    
    }
    
    
  }
  
  resultsTable<-do.call(cbind,resultsTable)
  
  #extract the ensembl gene_ids
  if (species=="Human"){
    endf <- GenomicFeatures::genes(EnsDb.Hsapiens.v79, return.type="DataFrame")
  } else if (species=="Rat"){
    endf <- GenomicFeatures::genes(EnsDb.Rnorvegicus.v79, return.type="DataFrame")
  } else if (species=="Cow"){
    endf <- GenomicFeatures::genes(EnsDb.Btaurus.v79, return.type="DataFrame")
  } else if (species=="Horse"){
    endf <- GenomicFeatures::genes(EnsDb.Ecaballus.v79, return.type="DataFrame")
  } else if (species=="Zebrafish"){
    endf <- GenomicFeatures::genes(EnsDb.Drerio.v79, return.type="DataFrame")
  } else if (species=="Pig"){
  endf <- GenomicFeatures::genes(EnsDb.Sscrofa.v79, return.type="DataFrame")
  } else {
    endf <- GenomicFeatures::genes(EnsDb.Mmusculus.v79, return.type="DataFrame")
  }
  
  #make a dataframe of gene_ids and names
  en2gene <- as.data.frame(endf[,c("gene_id","gene_name")])
  
  #merge the ids
  resultsTable<-merge(resultsTable,en2gene,by.x="row.names",by.y="gene_id")
  
  #tidy the table
  resultsTable<-resultsTable[,c(ncol(resultsTable),3:ncol(resultsTable)-1,1)]
  colnames(resultsTable)[ncol(resultsTable)]<-"ID"
  
  write.table(resultsTable,file=foldChangeTable,col.names=T,row.names=F,sep="\t",quote=F)


}
