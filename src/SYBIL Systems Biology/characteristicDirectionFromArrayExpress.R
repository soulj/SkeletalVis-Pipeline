#
# Uses the bioconductor library to perform differntial expression analysis
#
# 
# A expressionset/ affybatch object is expected
# 
# 
suppressPackageStartupMessages(library(RGalaxy))

characteristicDirectionFromArrayExpress <- function(
  inputfile = GalaxyInputFile(required=TRUE,formatFilter = "rdata"), 
  comparisonsTable=GalaxyInputFile(required=TRUE,formatFilter = "tabular"),
  platform=GalaxySelectParam(c("Affy","Affy-ST","1C-Agilent","2C-Agilent","Illumina","Illumina_Exp"),required=T), 
  annotationFile=GalaxySelectParam(c("bovine.db","hgu133a.db","hgu133plus2.db","hgug4112a.db","HsAgilentDesign026652.db","HsAgilentDesign026652.db",
                                     "hugene10sttranscriptcluster.db","hugene20sttranscriptcluster.db","lumiHumanIDMapping","lumiMouseIDMapping",
                                     "lumiRatIDMapping","mgug4122a.db","MmAgilentDesign026655.db", "moe430a.db","mogene10sttranscriptcluster.db","mogene20sttranscriptcluster.db",
                                     "mouse4302.db", "mouse430a2.db","porcine.db","ragene10sttranscriptcluster.db","ragene20sttranscriptcluster.db","rat2302.db","xlaevis.db","hgu95av2.db",
                                     "hgfocus.db","mgu74av2.db","u133x3p.db","primeview.db","AgilentMouse014868.db","AgilentMouse028005.db","AgilentRat014879.db","AgilentRat028279.db","AgilentHuman039494.db","ArrayXHuman.db","AgilentMouse074809.db","mta10transcriptcluster.db","htmg430pm.db","huex10sttranscriptcluster.db",
                                     "mgu74a.db","hugene21sttranscriptcluster.db","hta20transcriptcluster.db","AgilentMouse079303.db"),required=T),
  foldChangeOnly=GalaxyLogicalParam(),offset=GalaxyIntegerParam(0L),
  chrDirTable=GalaxyOutput("chrDirTable","tabular")) {
  
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(oligo))
  suppressPackageStartupMessages(library(affy))
  suppressPackageStartupMessages(library(gcrma))
  suppressPackageStartupMessages(library(affyQCReport))
  suppressPackageStartupMessages(library(affyPLM))
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(vsn))
  suppressPackageStartupMessages(library(sva))
  suppressPackageStartupMessages(library(annotationFile,character.only=TRUE))
  suppressPackageStartupMessages(library(GeoDE))
  
  regressSVs = function(y, mod, svaobj,  P=ncol(mod)) {
    X<-cbind(mod,svaobj$sv)
    Hat<-solve(t(X)%*%X)%*%t(X)
    beta<-(Hat%*%t(y))
    cleany<-y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
    return(cleany)
  }
  
  
  getResultsDataFrame<-function(comparisonNum,exprs,eset,comparisonsTable) {
    
    contrast <- comparisonsTable[comparisonNum,1]
    numerator <- comparisonsTable[comparisonNum,2]
    denominator <- comparisonsTable[comparisonNum,3]
    
    keep <- which(pData(eset)[[contrast]] %in% c(numerator,denominator))
    exprs <- exprs[,keep]
    eset <- eset[,keep]
    classes <- as.factor(ifelse(pData(eset)[[contrast]] %in% numerator,1,2))
    
    exprs <- as.data.frame(exprs)
    exprs <- cbind(rownames(exprs),exprs)
    
    keep <- na.omit(order(apply(exprs[,-1],1,var),decreasing = T)[1:20000])
    dim(exprs)
    exprs <- exprs[keep,]
    dim(exprs)
    chdir <- chdirAnalysis(exprs,classes,1,CalculateSig=F)
    # results <- data.frame(ID=rownames(chdir$chdirprops$chdir[[1]]),chrDir=chdir$chdirprops$chdir[[1]][,1])
    # results$sig <-ifelse(results$ID %in% names(chdir$results[[1]]),1,0)
    # results <- stack(chdir$results[[1]])[,2:1]
    # results$comparison <- paste(numerator,denominator,sep="vs")
    # results$comparisonNumber <- comparisonNum
    # colnames(results)[1:2] <- c("ID","chrDir")
    
    results.top500 <- stack(chdir$chdirprops$chdir[[1]][,1])[,2:1]
    results.top500 <- results.top500[ order(abs(results.top500[,2]),decreasing = T),]
    results.top500 <- results.top500[1:500,]
    results.top500$comparison <- paste(numerator,denominator,sep="vs")
    results.top500$comparisonNumber <- comparisonNum
    colnames(results.top500)[1:2] <- c("ID","chrDir")
    s
    return(results.top500)
    
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
  
  #### Preprocessing ####
  
  
  #get the object whatever it was named
  inputfile<-load(inputfile)
  inputfile <- get(inputfile)
  
  comparisonsTable<-read.delim(comparisonsTable,header=T,stringsAsFactors = F)
  
  if (foldChangeOnly == TRUE) {
    return("need replicates to perform characteristic direction analysis")
  }
  
  # if ("Paired2" %in% colnames(comparisonsTable)){
  #   
  #   return("characteristic direction cannot take into account paired samples")
  #   
  # }
  
  
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
    
    #don't aggregate the probe level to gene level since illumina contains one probe per gene
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
    
    rownames(exprs) <- rownames(eset$M)
    
    eset <- ExpressionSet(assayData = as.matrix(exprs),
                          phenoData = new("AnnotatedDataFrame", as.data.frame(targets)))
  
  
  }
  
  ###characteristic direction analysis#####
  
    #test for multiple factors to combine
  
  if (platform == "2C-Agilent"){
    
    comparisonsTable <- comparisonsTable[,c(5,4,3)]
    
    #make the comparisons
    resultsTable<- lapply(1:nrow(comparisonsTable), function(x) getResultsDataFrame(x,exprs,eset,comparisonsTable))
    resultsTable<-do.call(rbind,resultsTable)
    write.table(resultsTable,file=chrDirTable,col.names=T,row.names=F,sep="\t",quote=F)
    return()
    
  } else if ("Factor2" %in% colnames(comparisonsTable)){
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
      #make the design matrix
      designFormula <- as.formula(paste("~0 +", paste(factors, collapse=" + ")))
      design<-model.matrix(designFormula,data = factorValues)
      colnames(design)<-make.names(colnames(design))
      
    } else if ("Paired2" %in% colnames(comparisonsTable)){
      Factor1<-unique(as.character(comparisonsTable[,1]))
      Factor2<-unique(as.character(comparisonsTable[,2]))
      factors<-c(Factor1,Factor2)
      factors<-gsub("/| ",".",factors)
      colnames(pData(eset))<-gsub("/| ,:",".",colnames(pData(eset)))
      factorValues<-pData(eset)[factors]
      factorValues<-as.data.frame(apply(factorValues,2,as.factor))
      comparisonsTable<-comparisonsTable[,-1]
      comparisonsTable[,2]<-gsub("/| ,:",".",comparisonsTable[,2])
      #make the design matrix
      designFormula <- as.formula(paste("~0+", paste(factors[2:1], collapse=" + ")))
      design<-model.matrix(designFormula,data = factorValues)
      colnames(design)<-make.names(colnames(design))
      comparisonsTable$Paired2<-make.names(comparisonsTable$Paired2)
      pData(eset)[,Factor2]<-make.names(pData(eset)[,Factor2])
      
    } else {
      #define the experimental conditions (factors)
      factors<-as.character(unique(comparisonsTable[,1]))
      print(factors)
      factors<-gsub("/| ",".",factors)
      colnames(pData(eset))<-gsub("/| ,:",".",colnames(pData(eset)))
      comparisonsTable[,1]<-gsub("/| ,:",".",comparisonsTable[,1])
      factorValues<-pData(eset)[factors]
      if(any(grepl("+/+|-/-|+", factorValues[,1]))){
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "+/+",replacement = "WT",fixed = T)
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "-/-",replacement = "KO",fixed=T)
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "+/-",replacement = "HET",fixed=T)
        factorValues[,1]<-gsub(x = factorValues[,1],pattern = "+",replacement = ".",fixed=T)
      }
      if(any(grepl("^[[:digit:]]", factorValues[,1]))){
        factorValues<-as.data.frame(apply(factorValues,2,make.names))
      }
      factorValues<-as.data.frame(apply(factorValues,2,as.factor))
      #make the design matrix
      designFormula <- as.formula(paste("~0 +", paste(factors, collapse=" + ")))
      design<-model.matrix(designFormula,data = factorValues)
      colnames(design)<-make.names(colnames(design))
      
    }
  
  if (!"Paired2" %in% colnames(comparisonsTable)){
    
    #sva batch effect correction
    mod0 = model.matrix(~1,data=factorValues)
    svafit <- sva(exprs,mod = design,mod0=mod0)
    
    #regress out the batch effect
    if (!is.null(dim(svafit$sv))){
      exprs<-regressSVs(exprs,design,svafit)
    }
    
    if(any(grepl("+/+|-/-|+|-", pData(eset)[[factors]]))){
      pData(eset)[[factors]]<-gsub(x = pData(eset)[[factors]],pattern = "+/+",replacement = "WT",fixed = T)
      pData(eset)[[factors]]<-gsub(x = pData(eset)[[factors]],pattern = "-/-",replacement = "KO",fixed=T)
      pData(eset)[[factors]]<-gsub(x = pData(eset)[[factors]],pattern = "+/-",replacement = "HET",fixed=T)
      pData(eset)[[factors]]<-gsub(x = factorValues[,1],pattern = "+",replacement = ".",fixed=T)
      pData(eset)[[factors]]<-gsub(x = factorValues[,1],pattern = "-",replacement = ".",fixed=T)
      comparisonsTable[,2]<-gsub(x=comparisonsTable[,2],pattern = "-",replacement = ".",fixed=T)
      comparisonsTable[,3]<-gsub(x=comparisonsTable[,3],pattern = "-",replacement = ".",fixed=T)
      comparisonsTable[,2]<-gsub(x=comparisonsTable[,2],pattern = "+",replacement = ".",fixed=T)
      comparisonsTable[,3]<-gsub(x=comparisonsTable[,3],pattern = "+",replacement = ".",fixed=T)
      
    }
    if(any(grepl("^[[:digit:]]", pData(eset)[[factors]]))){
      
      pData(eset)[[factors]]<-make.names(pData(eset)[[factors]])
      
    }
  }

  #make the comparisons
  resultsTable<- lapply(1:nrow(comparisonsTable), function(x) getResultsDataFrame(x,exprs,eset,comparisonsTable))
  resultsTable<-do.call(rbind,resultsTable)
  write.table(resultsTable,file=chrDirTable,col.names=T,row.names=F,sep="\t",quote=F)
  
  
}

