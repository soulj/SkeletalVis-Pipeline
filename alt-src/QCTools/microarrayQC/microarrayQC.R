#!/usr/bin/env Rscript

## begin warning handler
withCallingHandlers({
  
  library(methods) # Because Rscript does not always do this
  
  options('useFancyQuotes' = FALSE)
  
  suppressPackageStartupMessages(library("optparse"))
  suppressPackageStartupMessages(library("RGalaxy"))
  
  
  option_list <- list()
  
  option_list$inputfile <- make_option('--inputfile', type='character')
  option_list$comparisonsTable <- make_option('--comparisonsTable', type='character')
  option_list$platform <- make_option('--platform', type='character')
  option_list$accessionNumber <- make_option('--accessionNumber', type='character')
  option_list$outputFile <- make_option('--outputFile', type='character')
  option_list$outputDirectory <- make_option('--outputDirectory', type='character')
  option_list$offset <- make_option('--offset', type='integer')


  opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("RGalaxy"))

microarrayQC <- function(inputfile,comparisonsTable,platform,accessionNumber,outputFile,outputDirectory,offset=0L) {

  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(oligo))
  suppressPackageStartupMessages(library(affy))
  suppressPackageStartupMessages(library(gcrma))
  suppressPackageStartupMessages(library(affyQCReport))
  suppressPackageStartupMessages(library(affyPLM))
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(vsn))
  suppressPackageStartupMessages(library(arrayQualityMetrics))
  #### Preprocessing ####
  
  if(platform=="2C-Agilent"){
    stop("2C-Agilent QC not possible at present")
  }
  
  
  #get the object whatever it was named
  inputfile<-load(inputfile)
  inputfile <- get(inputfile)
  
  comparisonsTable<-read.delim(comparisonsTable,header=T,stringsAsFactors = F)
  
  
  if (platform=="Affy"){
    
    ## normalize the data
    eset<-affy::rma(inputfile)


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

    
  } else if (platform == "Illumina"){
    
    ## normalize the data
    exprs(inputfile)<-exprs(inputfile)+offset
    eset<-normalize(inputfile,transfn="log")
    exprs(eset)<-log2(exprs(eset))

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
    
    #convert the normalised data back into an RG list for the PCA
    eset<-RG.MA(eset)
    expMat<-cbind(eset$R,eset$G)
    
  }
  
  
  arrayQualityMetrics(eset,outdir = outputDirectory,reporttitle = accessionNumber,force = T,intgroup = comparisonsTable[1,1])
  
  #rename the index.html to get it to match the galaxy generated outputFile name
  file.rename(paste0(outputDirectory,"/index.html"),outputFile)
  
  
}


params <- list()
for(param in names(opt))
{
  if (!param == "help")
    params[param] <- opt[param]
}

setClass("GalaxyRemoteError", contains="character")
wrappedFunction <- function(f)
{
  tryCatch(do.call(f, params),
           error=function(e) new("GalaxyRemoteError", conditionMessage(e)))
}


suppressPackageStartupMessages(library(RGalaxy))
do.call(microarrayQC, params)

## end warning handler
}, warning = function(w) {
  cat(paste("Warning:", conditionMessage(w), "\n"))
  invokeRestart("muffleWarning")
})
