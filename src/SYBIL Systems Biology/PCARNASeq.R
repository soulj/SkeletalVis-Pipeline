suppressPackageStartupMessages(library("RGalaxy"))


PCARNASeq<-function(txiData=GalaxyInputFile(required=T,formatFilter="rdata"),sampleTable=GalaxyInputFile(required=T,formatFilter="tabular"),species=GalaxySelectParam(c("Human","Mouse","Rat","Pig","Cow","Horse","Zebrafish"),required=T),technicalReplicates=GalaxyLogicalParam(),batchCorrect=GalaxyLogicalParam(),supervised=GalaxyLogicalParam(),pcaPlot=GalaxyOutput("PCA", "png"),geneInfluence=GalaxyOutput("geneInfluence","tabular")){
  
  #modification of DESeq2 plotPCA function to improve the figure for two experimental factors.
  PCA<- function(object, intgroup, ntop=30000,correctedMatrix)
  {
    if (is.null(correctedMatrix)){
      
      # calculate the variance for each gene
      rv <- rowVars(assay(object))
      
      # select the ntop genes by variance
      select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
      
      # perform a PCA on the data in assay(x) for the selected genes
      pca <- prcomp(t(assay(object)[select,]))
    
    } else {
     
      # perform a PCA on the data in assay(x) for the selected genes
      pca <- prcomp(t(correctedMatrix))
      
    }
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
    
    # assemble the data for the plot
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], intgroup.df, name=colnames(object))
    
    print(intgroup)
    print(length(intgroup))
    print(dim(intgroup.df))
    if (length(intgroup)>1){
      g<-ggplot(data=d, aes_string(x="PC1", y="PC2", color=intgroup.df[,1],shape=intgroup.df[,2])) + geom_point(size=3) + xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +       coord_fixed() + scale_colour_discrete(name=colnames(intgroup.df)[1])+scale_shape_discrete(name=colnames(intgroup.df)[2])
    } else {
      g<-ggplot(data=d, aes_string(x="PC1", y="PC2", color=intgroup.df[,1])) + geom_point(size=3) + xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +  coord_fixed() + scale_colour_discrete(name=colnames(intgroup.df)[1])
    }
    
    g<-g+cowplot::theme_cowplot()
    loadings<-as.data.frame(pca$rotation[,1:2])
    loadings$ID<-rownames(loadings)
    
    return (list(plot=g,loadings=loadings))
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
  
  load(txiData)
  
  sampleTable<-read.delim(sampleTable)
  
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
    
    pca <- PCA(varianceStabilizingTransformation(dds),intgroup=factors,correctedMatrix = correctedMatrix)
    
  } else {
    pca <- PCA(varianceStabilizingTransformation(dds),intgroup=factors,correctedMatrix = NULL)
  }
  
  
  ggsave(filename = pcaPlot,plot = pca$plot,device="png")
  
  #make a table of genes with the contribution to each of PCs
  loadings<-pca$loadings
  
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
  loadings<-merge(loadings,en2gene,by.x="ID",by.y="gene_id")
  
  #tidy the table
  loadings<-loadings[,c(1,4,2,3)]
  loadings<-loadings[ order(abs(loadings$PC1),decreasing = T),]
  
  write.table(loadings,file=geneInfluence,col.names=T,row.names=F,sep="\t",quote=F)
  
}
