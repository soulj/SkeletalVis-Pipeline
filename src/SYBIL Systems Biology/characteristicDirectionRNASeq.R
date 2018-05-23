suppressPackageStartupMessages(library("RGalaxy"))


characteristicDirectionRNASeq<-function(txiData=GalaxyInputFile(required=T,formatFilter="rdata"),
sampleTable=GalaxyInputFile(required=T,formatFilter="tabular"),
comparisonsTable=GalaxyInputFile(required=T,formatFilter="tabular"),
species=GalaxySelectParam(c("Human","Mouse","Rat","Pig","Cow","Horse","Zebrafish"),required=T),
foldChangeOnly=GalaxyLogicalParam(),
technicalReplicates=GalaxyLogicalParam(),batchCorrect=GalaxyLogicalParam(),
supervised=GalaxyLogicalParam(),chrDirTable=GalaxyOutput("chrDirTable","tabular")){

    getResultsDataFrame<-function(dds,comparisonNum,comparisonsTable,correctedMatrix) {
      
      #define the comparison to be made
      contrast <- comparisonsTable[comparisonNum,1]
      numerator <- comparisonsTable[comparisonNum,2]
      denominator <- comparisonsTable[comparisonNum,3]
      
      #subset the counts matrix to keep just the samples in the comparisons
      keep <- which(colData(dds)[[contrast]] %in% c(numerator,denominator))
      
      if (is.null(correctedMatrix)){
        
        exprs <- as.matrix(assay(varianceStabilizingTransformation(dds)))
        
      } else {

       exprs <- correctedMatrix
        
      }
      
      exprs <- exprs[,keep]
      dds <- dds[,keep]
      classes <- as.factor(ifelse(colData(dds)[[contrast]] %in% numerator,1,2))
      
      exprs <- as.data.frame(exprs)
      exprs <- cbind(rownames(exprs),exprs)
      
      keep <- na.omit(order(apply(exprs[,-1],1,var),decreasing = T)[1:20000])
      dim(exprs)
      exprs <- exprs[keep,]
      dim(exprs)
      
      
      chdir <- chdirAnalysis(exprs,classes,1,CalculateSig=F)
      
      results.top500 <- stack(chdir$chdirprops$chdir[[1]][,1])[,2:1]
      results.top500 <- results.top500[ order(abs(results.top500[,2]),decreasing = T),]
      results.top500 <- results.top500[1:500,]
      results.top500$comparison <- paste(numerator,denominator,sep="vs")
      results.top500$comparisonNumber <- comparisonNum
      colnames(results.top500)[1:2] <- c("ID","chrDir")
      
      return(results.top500)
    }

    
  regressSVs = function(y, mod, svaobj,  P=ncol(mod)) {
    X<-cbind(mod,svaobj$sv)
    Hat<-solve(t(X)%*%X)%*%t(X)
    beta<-(Hat%*%t(y))
    cleany<-y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
    return(cleany)
  }
  
  suppressPackageStartupMessages(library("DESeq2"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("cowplot"))
  suppressPackageStartupMessages(library("genefilter"))
  suppressPackageStartupMessages(library("EnsDb.Hsapiens.v79"))
  suppressPackageStartupMessages(library("EnsDb.Mmusculus.v79"))
  suppressPackageStartupMessages(library("EnsDb.Rnorvegicus.v79"))
  suppressPackageStartupMessages(library("EnsDb.Ecaballus.v79"))
  suppressPackageStartupMessages(library("EnsDb.Drerio.v79"))
  suppressPackageStartupMessages(library("EnsDb.Btaurus.v79"))
  suppressPackageStartupMessages(library("EnsDb.Sscrofa.v79"))
  suppressPackageStartupMessages(library("sva"))
  suppressPackageStartupMessages(library("GeoDE"))
  
  if (foldChangeOnly == TRUE) {
    return("need replicates to perfrom characteristic direction analysis")
  }
  
  load(txiData)
  
  sampleTable<-read.delim(sampleTable)
  comparisonsTable <- read.delim(comparisonsTable,stringsAsFactors = F)
  
  
  #remove the file info
  sampleTable<-sampleTable[,colnames(sampleTable)!="File"]
  
  #match the txi dataframes to the order of the sampleTable sample names
  txi$counts <- txi$counts[ ,match(as.character(sampleTable[,1]),colnames(txi$counts))]
  txi$abundance <- txi$abundance[ ,match(as.character(sampleTable[,1]),colnames(txi$abundance))]
  txi$length <- txi$length[ ,match(as.character(sampleTable[,1]),colnames(txi$length))]
  
  #get the design formula from the columns in the sampleTable (using the same order)
  factors<-colnames(sampleTable)[-1]
  designFormula <- as.formula(paste("~", paste(factors, collapse=" + ")))
  
  #perform normalisation and create the PCA plot
  #if technical replicates need collapsing
  if (technicalReplicates){
    dds<-DESeqDataSetFromTximport(txi = txi,colData = sampleTable,design = designFormula,ignoreRank = T)
    dds<-collapseReplicates(dds,paste0(sampleTable[,factors[1]],sampleTable[,factors[2]]))
    design(dds)<-as.formula(paste("~", paste(factors[-1], collapse=" + ")))
  } else {
    dds<-DESeqDataSetFromTximport(txi = txi,colData = sampleTable,design = designFormula)
  }
  
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
    
    correctedMatrix <- regressSVs(dat,mod,svseq)
    correctedMatrix <- correctedMatrix + abs(min(correctedMatrix))
    
    dds<-DESeqDataSetFromMatrix(round(correctedMatrix),colData = sampleTable,design = designFormula)
    correctedMatrix<-as.matrix(assay(varianceStabilizingTransformation(dds)))
    
    resultsTable <-  lapply(1:nrow(comparisonsTable),function(x) getResultsDataFrame(dds,comparisonNum = x,comparisonsTable = comparisonsTable, correctedMatrix = correctedMatrix))
  
  } else {
    resultsTable <-  lapply(1:nrow(comparisonsTable),function(x) getResultsDataFrame(dds,comparisonNum = x,comparisonsTable = comparisonsTable, correctedMatrix = NULL))
    
  }
  
  resultsTable<-do.call(rbind,resultsTable)

  
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
  resultsTable<-merge(resultsTable,en2gene,by.x="ID",by.y="gene_id")
  resultsTable <- resultsTable[ order(resultsTable$comparisonNumber,rev(abs(resultsTable$chrDir))),]

  
  write.table(resultsTable,file=chrDirTable,col.names=T,row.names=F,sep="\t",quote=F)
  
}
