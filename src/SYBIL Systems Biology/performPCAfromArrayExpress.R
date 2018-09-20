suppressPackageStartupMessages(library(RGalaxy))

performPCAfromArrayExpress <- function(
  inputfile = GalaxyInputFile(required=TRUE,formatFilter = "rdata"), 
  comparisonsTable=GalaxyInputFile(required=TRUE,formatFilter = "tabular"),
  platform=GalaxySelectParam(c("Affy","Affy-ST","1C-Agilent","2C-Agilent","Illumina"),required=T), 
  annotationFile=GalaxySelectParam(c("bovine.db","hgu133a.db","hgu133plus2.db","hgug4112a.db","HsAgilentDesign026652.db","HsAgilentDesign026652.db",
                                     "hugene10sttranscriptcluster.db","hugene20sttranscriptcluster.db","lumiHumanIDMapping","lumiMouseIDMapping",
                                     "lumiRatIDMapping","mgug4122a.db","MmAgilentDesign026655.db", "moe430a.db","mogene10sttranscriptcluster.db","mogene20sttranscriptcluster.db",
                                     "mouse4302.db", "mouse430a2.db","porcine.db","ragene10sttranscriptcluster.db","ragene20sttranscriptcluster.db","rat2302.db","xlaevis.db","hgu95av2.db",
                                     "hgfocus.db","mgu74av2.db","u133x3p.db","primeview.db","AgilentMouse014868.db","AgilentMouse028005.db","AgilentRat014879.db","AgilentRat028279.db","AgilentHuman039494.db","ArrayXHuman.db","AgilentMouse074809.db","mta10transcriptcluster.db","htmg430pm.db","huex10sttranscriptcluster.db",
                                     "mgu74a.db","hugene21sttranscriptcluster.db","hta20transcriptcluster.db","AgilentMouse079303.db"),required=T),
  foldChangeOnly=GalaxyLogicalParam(),offset=GalaxyIntegerParam(0L) , pcaPlot = GalaxyOutput("PCAplot", "png"),geneInfluence=GalaxyOutput("geneInfluence", "tabular")) {
  
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(oligo))
  suppressPackageStartupMessages(library(affy))
  suppressPackageStartupMessages(library(gcrma))
  suppressPackageStartupMessages(library(affyQCReport))
  suppressPackageStartupMessages(library(affyPLM))
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(vsn))
  suppressPackageStartupMessages(library(sva))
  suppressPackageStartupMessages(library(cowplot))
  suppressPackageStartupMessages(library(annotationFile,character.only=TRUE))
  
  
  #calculate PCA from a expressionSet type object
  PCA.eset<- function(exprs,eset, expFactors, ntop=30000)
  {
    # calculate the variance for each gene
    rv <- rowVars(exprs)
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    # perform a PCA on the data for the selected genes
    pca <- prcomp(t(exprs[select,]))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    
    # assemble the data for the plot
    intgroup.df <- as.data.frame(pData(eset)[,expFactors], drop=FALSE)
    
    #fix the names in case there are spaces
    intgroup.df<-as.data.frame(apply(intgroup.df,2,function(x) paste("`", x, "`", sep="")))
    colnames(intgroup.df)<-colnames(pData(eset)[expFactors])
    
    
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],intgroup.df )
    
    
    if (length(expFactors)>1){
      g<-ggplot(data=d, aes_string(x="PC1", y="PC2", color=intgroup.df[,1],shape=intgroup.df[,2])) + geom_point(size=3) + xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +       coord_fixed() + scale_colour_discrete(name=colnames(intgroup.df)[1])+scale_shape_discrete(name=colnames(intgroup.df)[2])
      
    } else {
      g<-ggplot(data=d, aes_string(x="PC1", y="PC2", color=intgroup.df[,1])) + geom_point(size=3) + xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +  coord_fixed() + scale_colour_discrete(name=colnames(intgroup.df)[1])
      }
    #g<-g+cowplot::theme_cowplot()+theme(legend.position="bottom")+guides(col = guide_legend(nrow = length(unique(intgroup.df[,1]))))
    g<-g+cowplot::theme_cowplot()
    loadings<-as.data.frame(pca$rotation[,1:2])
    loadings$ID<-rownames(loadings)
    
    loadings<-loadings[,c(3,1,2)]
    loadings<-loadings[ order(abs(loadings$PC1),decreasing = T),]
    
    
    return (list(plot=g,loadings=loadings))
  }
  
  
  #calculate PCA from a expressionSet type object
  PCA.matrix<- function(expMatrix,phenoData,expFactors, ntop=30000)
  {
    # calculate the variance for each gene
    rv <- rowVars(expMatrix)
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    # perform a PCA on the data for the selected genes
    pca <- prcomp(t(expMatrix[select,]))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    
    # assemble the data for the plot
    intgroup.df <- as.data.frame(phenotypes,drop=F)
    
    #fix the names in case there are spaces
    intgroup.df<-as.data.frame(apply(intgroup.df,1,function(x) paste("`", x, "`", sep="")))
    colnames(intgroup.df)<-colnames(phenotypes)
    
    
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],intgroup.df )
    
    
    if (length(expFactors)>1){
      g<-ggplot(data=d, aes_string(x="PC1", y="PC2", color=intgroup.df[,1],shape=intgroup.df[,2])) + geom_point(size=3) + xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +       coord_fixed() + scale_colour_discrete(name=colnames(intgroup.df)[1])+scale_shape_discrete(name=colnames(intgroup.df)[2])
    } else {
      g<-ggplot(data=d, aes_string(x="PC1", y="PC2", color=intgroup.df[,1])) + geom_point(size=3) + xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +  coord_fixed() + scale_colour_discrete(name=colnames(intgroup.df)[1])
    }
    g<-g+cowplot::theme_cowplot()
    
    loadings<-as.data.frame(pca$rotation[,1:2])
    loadings$ID<-rownames(loadings)
    
    loadings<-loadings[,c(3,1,2)]
    loadings<-loadings[ order(abs(loadings$PC1),decreasing = T),]
    
    
    return (list(plot=g,loadings=loadings))
  }
  
  
  annotateProbes<-function(exprs,annotationFile) {
    
    #annotate the probes
    annotationFile<-gsub(".db",replacement = "",annotationFile)
    geneIDs <- na.omit(stack(mget(as.character(rownames(exprs)), get(paste(annotationFile,"SYMBOL",sep="")), ifnotfound = NA)))
    exprs<-merge(exprs,geneIDs,by.x="row.names",by.y="ind")
    exprs$Row.names=NULL
    
    #aggregate the probe level to gene level
    exprs<-aggregate(exprs[,-ncol(exprs)],by=list(exprs$values),FUN=median,na.rm=TRUE)
    rownames(exprs)<-exprs$Group.1
    exprs[,1]<-NULL
    exprs<-as.matrix(exprs)
    return(exprs)
  }
  
  regressSVs = function(y, mod, svaobj,  P=ncol(mod)) {
    X<-cbind(mod,svaobj$sv)
    Hat<-solve(t(X)%*%X)%*%t(X)
    beta<-(Hat%*%t(y))
    cleany<-y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
    return(cleany)
  }
  
  #### Preprocessing ####
  
  # if(platform=="2C-Agilent"){
  #   stop("2C-Agilent PCA not possible at present")
  # }
  
  
  #get the object whatever it was named
  inputfile<-load(inputfile)
  inputfile <- get(inputfile)
  
  comparisonsTable<-read.delim(comparisonsTable,header=T,stringsAsFactors = F)
  
  
  if (platform=="Affy"){
    
    ## normalize the data
    eset<-affy::rma(inputfile)
    exprs<-exprs(eset)
    exprs<-annotateProbes(exprs,annotationFile)
    
    
  } else if (platform=="Affy-ST"){
    
    suppressPackageStartupMessages(library(pd.hugene.1.0.st.v1))
    suppressPackageStartupMessages(library(pd.hugene.2.0.st))
    suppressPackageStartupMessages(library(pd.mogene.1.0.st.v1))
    suppressPackageStartupMessages(library(pd.mogene.2.0.st))
    suppressPackageStartupMessages(library(pd.ragene.1.0.st.v1))
    suppressPackageStartupMessages(library(pd.ragene.2.0.st))
    suppressPackageStartupMessages(library(pd.mta.1.0))
    suppressPackageStartupMessages(library(pd.huex.1.0.st.v2))
    suppressPackageStartupMessages(library(pd.hugene.2.1.st))
    suppressPackageStartupMessages(library(pd.hta.2.0))
    
    ## normalize the data
    eset<-oligo::rma(inputfile)
    exprs<-exprs(eset)
    exprs<-annotateProbes(exprs,annotationFile)
    
    
  }else if (platform == "Illumina"){
    
    ## normalize the data
    exprs(inputfile)<-exprs(inputfile)+offset
    eset<-normalize(inputfile,transfn="log")
    exprs<-log2(exprs(eset))
    
    #annotate the probes by the nuID method
    if (annotationFile=="lumiMouseIDMapping"){
      rownames(exprs)<-as.data.frame(IlluminaID2nuID(rownames(exprs),lib.mapping = 'lumiMouseIDMapping',species="Mouse"))$Symbol
    } else{
      rownames(exprs)<-as.data.frame(IlluminaID2nuID(rownames(exprs),lib.mapping = 'lumiHumanIDMapping',species="Human"))$Symbol
    }
    
    #aggregate the probe level to gene level
    exprs<-aggregate(exprs,by=list(rownames(exprs)),FUN=median,na.rm=TRUE)
    rownames(exprs)<-exprs$Group.1
    exprs[,1]<-NULL
    exprs<-as.matrix(exprs)
    
    
  }   else if (platform == "1C-Agilent"){
    
    ## normalize
    eset<-limma::normalizeBetweenArrays(inputfile)
    
    #aggregate the duplicated probes by the median intensity
    rownames(eset$E)<-eset$genes$ProbeName
    rownames(eset$targets)<-colnames(eset$E)
    eset$E <-aggregate(eset$E,by=list(eset$genes$ProbeName),FUN=median,na.rm=TRUE)
    rownames(eset$E)<-eset$E$Group.1
    eset$E$Group.1=NULL
    
    #convert to eset to work with general code
    eset <- ExpressionSet(assayData = as.matrix(eset$E),
                          phenoData = new("AnnotatedDataFrame", eset$targets))
    
    exprs<-exprs(eset)
    exprs<-annotateProbes(exprs,annotationFile)
    
  }   else if (platform == "2C-Agilent"){
    
      ## normalize within and between arrays
      eset <- limma::backgroundCorrect(inputfile, method="normexp", offset=50)
      eset<-limma::normalizeWithinArrays(eset,method="loess")
      eset<-limma::normalizeBetweenArrays(eset,method="Aquantile")
      
      #aggregate the duplicated gene probes
      annotationFile<-gsub(".db",replacement = "",annotationFile)
      geneIDs <- na.omit(stack(mget(as.character(eset$genes$ProbeName), get(paste(annotationFile,"SYMBOL",sep="")), ifnotfound = NA)))
      
      eset<-eset[ which(eset$genes$ProbeName %in% geneIDs$ind),]
      eset$genes$GeneName<-geneIDs$values
      eset<-avereps(eset, eset$genes$GeneName)
      
      Cy3<-comparisonsTable[1,"Cy3"]
      Cy5<-comparisonsTable[1,"Cy5"]
      eset$targets<-data.frame(FileName=colnames(eset$M),Cy5=gsub(' +',' ',eset$targets[,Cy5]),Cy3=gsub(' +',' ',eset$targets[,Cy3]))
      
      exprs <- exprs.MA(eset)
      targets <- data.frame(as.vector(t(cbind(as.character(eset$targets$Cy5),as.character(eset$targets$Cy3)))))
      colnames(targets)<-comparisonsTable$Factor[1]

      eset <- ExpressionSet(assayData = as.matrix(exprs),
                            phenoData = new("AnnotatedDataFrame", as.data.frame(targets)))
    }
  
  ##PCA
  
  #2 colour arrays require special treatment
  if (platform == "2C-Agilent"){
    
    pca<-PCA.eset(exprs,eset,comparisonsTable$Factor[1])

    
  } else if (foldChangeOnly == TRUE) {
    
    #test for multiple factors to combine
    if ("Factor2" %in% colnames(comparisonsTable)){
      Factor1<-unique(as.character(comparisonsTable[,1]))
      Factor2<-unique(as.character(comparisonsTable[,2]))
      factors<-paste(Factor1,Factor2,sep=".")
      pData(eset)[factors]<-paste(pData(eset)[,Factor1],pData(eset)[,Factor2],sep=".")
      factors<-gsub("/| ",".",factors)
      colnames(pData(eset))<-gsub("/| ",".",colnames(pData(eset)))
      factorValues<-pData(eset)[factors]
      factorValues<-as.data.frame(apply(factorValues,2,as.factor))
      comparisonsTable[,1]<-paste(comparisonsTable[,1],comparisonsTable[,2],sep=".")
      comparisonsTable[,1]<-gsub("/| ",".",comparisonsTable[,1])
      comparisonsTable<-comparisonsTable[,-2]
    }
    
    #make the fold change comparisons
    factors<-as.character(unique(comparisonsTable[,1]))
    pca<-PCA.eset(exprs,eset,factors)
    
  }  else {
    
    #test for multiple factors to combine
    if ("Factor2" %in% colnames(comparisonsTable)){
      Factor1<-unique(as.character(comparisonsTable[,1]))
      Factor2<-unique(as.character(comparisonsTable[,2]))
      factors<-paste(Factor1,Factor2,sep=".")
      pData(eset)[factors]<-paste(pData(eset)[,Factor1],pData(eset)[,Factor2],sep=".")
      factors<-gsub("/| ",".",factors)
      colnames(pData(eset))<-gsub("/| ",".",colnames(pData(eset)))
      factorValues<-pData(eset)[factors]
      if(any(grepl("+/+|-/-", factorValues[,1]))){
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "+/+",replacement = "WT",fixed = T)
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "-/-",replacement = "KO",fixed=T)
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "+/-",replacement = "HET",fixed=T)
      }
      factorValues<-as.data.frame(apply(factorValues,2,as.factor))
      comparisonsTable[,1]<-paste(comparisonsTable[,1],comparisonsTable[,2],sep=".")
      comparisonsTable[,1]<-gsub("/| ",".",comparisonsTable[,1])
      comparisonsTable<-comparisonsTable[,-2]
      
    } else if ("Paired1" %in% colnames(comparisonsTable)){
      Factor1<-unique(as.character(comparisonsTable[,1]))
      Factor2<-unique(as.character(comparisonsTable[,2]))
      factors<-c(Factor1,Factor2)
      factors<-gsub("/| ",".",factors)
      colnames(pData(eset))<-gsub("/| ",".",colnames(pData(eset)))
      factorValues<-pData(eset)[factors]
      factorValues<-as.data.frame(apply(factorValues,2,as.factor))
      comparisonsTable<-comparisonsTable[,-1]
      comparisonsTable[,2]<-gsub("/| ",".",comparisonsTable[,2])
      #make the design matrix
      designFormula <- as.formula(paste("~", paste(factors, collapse=" + ")))
      design<-model.matrix(designFormula,data = factorValues)
      colnames(design)<-make.names(colnames(design))
      
    } else {
      #define the experimental conditions (factors)
      factors<-as.character(unique(comparisonsTable[,1]))
      print(factors)
      factors<-gsub("/| ",".",factors)
      colnames(pData(eset))<-gsub("/| ",".",colnames(pData(eset)))
      comparisonsTable[,1]<-gsub("/| ",".",comparisonsTable[,1])
      factorValues<-pData(eset)[factors]
      if(any(grepl("^[[:digit:]]", factorValues[,1]))){
        factorValues<-as.data.frame(apply(factorValues,2,make.names))
      }
      if(any(grepl("+/+|-/-", factorValues[,1]))){
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "+/+",replacement = "WT",fixed = T)
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "-/-",replacement = "KO",fixed=T)
      }
      
      factorValues<-as.data.frame(apply(factorValues,2,as.factor))
      
    }
    
    #make the design matrix
    designFormula <- as.formula(paste("~0 +", paste(factors, collapse=" + ")))
    design<-model.matrix(designFormula,data = factorValues)
    colnames(design)<-make.names(colnames(design))
    
    if (!"Paired2" %in% colnames(comparisonsTable)){
    
      #add sva batch effect correction
      mod0 = model.matrix(~1,data=factorValues)
      svafit <- sva(exprs,mod = design,mod0=mod0)
      
      #regress the svs from the exprs matrix to make corrected matrix
      exprs<-regressSVs(exprs,design,svafit)
    }
    ## perform the pca
    pca<-PCA.eset(exprs,eset,factors)
    
  }
  
  #make a table of genes with the contribution to each of PCs
  loadings<-pca$loadings
  write.table(loadings,file=geneInfluence,col.names=T,row.names=F,sep="\t",quote=F)
  ggsave(filename = pcaPlot,plot = pca$plot,device="png",height = 7,width=9)
}
