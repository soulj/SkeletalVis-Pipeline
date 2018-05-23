suppressPackageStartupMessages(library("RGalaxy"))


#Function to download raw data from GEO given an accession number
getMicroarrayExpressionData<-function(accessionNumber=GalaxyCharacterParam(),
                                             platform=GalaxySelectParam(c("Affy","Affy-ST","1C-Agilent","2C-Agilent","Illumina")),
                                             numberFileRemove=GalaxyIntegerParam(0L), 
                                             grepExpression=GalaxyLogicalParam(),
                                             grepString=GalaxyCharacterParam(),
                                             numberLinesSkip=GalaxyIntegerParam(0L),
                                             split=GalaxyLogicalParam(checked=F),
                                             splitField=GalaxyCharacterParam(),
                                             splitSep=GalaxyCharacterParam(" "),
                                             splitPos=GalaxyIntegerParam(1L), 
                                             remove=GalaxyLogicalParam(checked=F),
                                             removeSample=GalaxyCharacterParam(),
                                             site=GalaxySelectParam(c("GEO","ArrayExpress"),required=T), 
                                             expressionData=GalaxyOutput("expressionData","rdata")) {
  
  suppressPackageStartupMessages(library("GEOquery"))
  suppressPackageStartupMessages(library("limma"))
  suppressPackageStartupMessages(library("oligo"))
  suppressPackageStartupMessages(library("simpleaffy"))
  suppressPackageStartupMessages(library("lumi"))
  
  #get the raw data files from arrayexpress`
  getArrayExpressSuppFiles<-function(accessionNumber){
    
    baseurl<-"https://www.ebi.ac.uk/arrayexpress/files/"
    endurl<-".raw.1.zip"
    downloadUrl<-paste(baseurl,accessionNumber,"/",accessionNumber,endurl,sep="")
    destFile<-paste(accessionNumber,".zip",sep="")
    
    # ebi can be temperamental, try to access 5 times before failing
    attempt <- 1
    while ((is.na(file.size(destFile)) || file.size(destFile) < 1000) && attempt <= 5) {
      if (attempt > 1) {
        Sys.sleep(60)
      }
      attempt <- attempt + 1
      try (
        download.file(downloadUrl,destFile)
      )
    }
    if (is.na(file.size(destFile)) || file.size(destFile) < 1000) {
      stop(sprintf("Failed to access raw data file for accession number %s", accessionNumber))
    }
    return(destFile)
  }
  #get the phenotype data file from arrayexpress
  getArrayExpresspData<-function(accessionNumber){
    
    baseurl<-"https://www.ebi.ac.uk/arrayexpress/files/"
    endurl<-".sdrf.txt"
    downloadUrl<-paste(baseurl,accessionNumber,"/",accessionNumber,endurl,sep="")
    destFile<-paste(accessionNumber,".txt",sep="")
    # ebi can be temperamental, try to access 5 times before failing
    attempt <- 1
    while ((is.na(file.size(destFile)) || file.size(destFile) < 1000) && attempt <= 5) {
      if (attempt > 1) {
        Sys.sleep(60)
      }
      attempt <- attempt + 1
      try (
        download.file(downloadUrl,destFile)
      )
    }
    if (is.na(file.size(destFile)) || file.size(destFile) < 1000) {
      stop(sprintf("Failed to access phenotype data for accession number %s", accessionNumber))
    }
    pdata<-read.delim(destFile,stringsAsFactors = F)
    
    return(pdata)
  }
  
  # Get phenotype data from GEO
  getGEOPhenotypeData<-function(accesionNumber){
    # GEO can be temperamental, try to access phenotype 5 times before failing
    eset <- try (getGEO(accessionNumber)[[1]])
    attempt <- 2
    
    while (inherits(eset, "try-error") && attempt <= 5) {
      Sys.sleep(60)
      attempt <- attempt + 1
      eset <- try (getGEO(accessionNumber)[[1]])
    }
    if (inherits(eset, "try-error")) {
      stop(sprintf("Failed to access GEO phenotype data for accession number %s", accessionNumber))
    }
    
    if (platform != "2C-Agilent"){
    columnNames <- sapply(strsplit(x = as.character(colnames(pData(eset))),split=":"),"[[",1)
    colnames(pData(eset)) <- make.names(columnNames,unique = T)
    pData(eset) <-  as.data.frame(apply(pData(eset),2,trimws))
    }
    
    return(eset)
  }
  
  
  #download the compressed archive of the raw data
  if (site=="GEO"){
    
      # GEO can be temperamental, try to access 5 times before failing
      filePaths <- try(getGEOSuppFiles(accessionNumber))
      
      attempt <- 2
      
      while (inherits(filePaths, "try-error") && attempt <= 5) {
        Sys.sleep(60)
        attempt <- attempt + 1
        filePaths <- try(getGEOSuppFiles(accessionNumber))
      }
      if (inherits(filePaths, "try-error")) {
        stop(sprintf("Failed to access GEO raw data for accession number %s", accessionNumber))
      }
      untar(tarfile = grep(pattern = ".tar",x = rownames(filePaths),value = T),exdir = accessionNumber)
      file.remove(grep(pattern = ".tar",x = rownames(filePaths),value = T))
  }  else {
    
    filePath<-getArrayExpressSuppFiles(accessionNumber)
    unzip(filePath,exdir = accessionNumber)
    file.remove(filePath)
    
  }
  
  # Transforms raw data into R format for analysis
  if (platform == "1C-Agilent") {
    
    if(site=="GEO"){
      
      #download the raw text files and removes all non raw data files
      files<-list.files(accessionNumber,pattern = ".txt")
      if (numberFileRemove>0){
        files<-files[-1:-numberFileRemove]
      }
      
      setwd(accessionNumber)
      
      #decompress
      sapply(files,gunzip)
      files<-list.files(pattern = ".txt")
      if (numberFileRemove>0){
        files<-files[-1:-numberFileRemove]
      }
      
      #read files with limma
      expData <- read.maimages(path = ".",files = files,source="agilent",green.only = T)
      
      eset<-getGEOPhenotypeData(accessionNumber)
      
      targets<-pData(eset)
      expData$targets<-targets
      
      setwd("..")
      
    } else{
      setwd(accessionNumber)
      files<-list.files(pattern = ".txt")
      
      if (numberFileRemove>0){
        files<-files[-1:-numberFileRemove]
      }
      
      
      #read files with limma
      expData <- read.maimages(path = ".",files = files,source="agilent",green.only = T)
      
      #get the phenotype data
      pdata<-getArrayExpresspData(accessionNumber)
      
      expData<-expData[,sort(colnames(expData$E))]
      pdata<-pdata[ order(as.character(pdata$Array.Data.File)),]
    
      expData$targets<-pdata
      
      setwd("..")
    }
    
  } else  if (platform == "2C-Agilent") {
    
    setwd(accessionNumber)
    files<-list.files(pattern = ".txt")
    
    if (numberFileRemove>0){
      files<-files[-1:-numberFileRemove]
    }
    
    
    #read files with limma
    expData <- read.maimages(path = ".",files = files,source="agilent",green.only = F)
    
    #get the phenotype data
    if(site=="GEO"){
      
      eset<-getGEOPhenotypeData(accessionNumber)
      targets<-pData(eset)
      expData$targets<-targets
      
    } else {
      pdata<-getArrayExpresspData(accessionNumber)
      expData$targets<-pdata
      
      n<-nrow(expData$targets)
      expData$targets<-merge(expData$targets[1:(n/2),],expData$targets[(n/2)+1:n,],by="Hybridization.Name")
      
      
    }
    
    setwd("..")
    
    
  } else  if (platform == "Affy-ST") {
    
    #get the CEL files
    files<-list.files(accessionNumber,pattern = ".CEL",ignore.case = T)
    
    setwd(accessionNumber)
    
    #read the data
    expData <- oligo::read.celfiles(files) 
    
    if(site=="GEO"){
    
      eset<-getGEOPhenotypeData(accessionNumber)
      #replace the pdata
      rownames(pData(eset))<-rownames(pData(expData))
      pData(expData)<-pData(eset)
    
    } else {
      
      pdata<-getArrayExpresspData(accessionNumber)
      
      expData<-expData[,sort(sampleNames(expData))]
      pdata<-pdata[ order(as.character(pdata$Array.Data.File)),]
      #replace the pdata
      rownames(pdata)<-rownames(pData(expData))
      pData(expData)<-pdata
      
    }
    
    
    setwd("..")
    
    
    
  }  else  if (platform == "Illumina") {
    
    if(site=="GEO"){
      
      #get the illumina files
      files<-list.files(accessionNumber,pattern = accessionNumber)
      
      setwd(accessionNumber)
      
      #decompress
      sapply(files,gunzip)
      file<-list.files(pattern = accessionNumber)
      if(length(file)>1){
        file<-file[length(file)]
      }
      
      #make the files into an eset if the data just contains expression values
      
      data<-read.delim(file,skip=numberLinesSkip)
      
      if (grepExpression==T){
        exprs<-data[,grep(grepString,colnames(data))] 
        
      } else {
        exprs<-data[,seq(from=2,to=ncol(data),by=2)]
        exprs <-  exprs[,colSums(is.na( exprs))<nrow( exprs)]
        
      }
      
      
      rownames(exprs)<-data[,1]
      
      eset<-getGEOPhenotypeData(accessionNumber)
      rownames(pData(eset))<-colnames(exprs)
      expData <- ExpressionSet(assayData = as.matrix(exprs),
                               phenoData = phenoData(eset))
      
      setwd("..")
    } else{
      
      file<-list.files(accessionNumber,pattern =".txt",full.names = T)
      expData<-lumiR(file)
      
      #get the pData
      pdata<-getArrayExpresspData(accessionNumber)
      
      #fix the row order
      expData<-expData[,sort(sampleNames(expData))]
      pdata<-pdata[ order(as.character(pdata$Array.Data.File)),]
      
      pData(expData)<-pdata
      rownames(pData(expData))<-colnames(exprs(expData))
      
      expData <- ExpressionSet(assayData = as.matrix(exprs(expData)),
                               phenoData = phenoData(expData))
    }
    
    
  }  else  if (platform == "Affy") {
    
    #get the CEL files
    files<-list.files(accessionNumber,pattern = ".CEL",ignore.case = T)
    
    setwd(accessionNumber)
    
    #read the data
    expData <- ReadAffy(filenames = files)
    
    if(site=="GEO"){
      eset<-getGEOPhenotypeData(accessionNumber)
      #replace the pdata
      rownames(pData(eset))<-rownames(pData(expData))
      pData(expData)<-pData(eset)
    }
    else {
      pdata<-getArrayExpresspData(accessionNumber)
      
      #fix the row order
      expData<-expData[,sort(sampleNames(expData))]
      
      pdata<-pdata[ order(as.character(pdata$Array.Data.File)),]

      #replace the pdata
      rownames(pdata)<-rownames(pData(expData))
      
      pData(expData)<-pdata
      
    }
    if(remove==TRUE){
      expData<-expData[,!sampleNames(expData)==removeSample]
    }
    
    setwd("..")
  }
  
  #tidy the pdata if needed
  if(split==TRUE){
    if (platform=="1C-Agilent"){
      expData$targets[,splitField]<-trimws(sapply(strsplit(as.character(expData$targets[,splitField]),split = splitSep),"[[",splitPos))
      
    } else{
      pData(expData)[splitField]<-trimws(sapply(strsplit(as.character(pData(expData)[,splitField]),split = splitSep),"[[",splitPos))
    }
    
  }
  
  
  save(expData,file = expressionData)
}




