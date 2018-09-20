library("NMF")
library("RColorBrewer")
library("preprocessCore")
library("dplyr")
library("org.Hs.eg.db")
library("feather")
library("taRifx")
library("purrr")

completed <- read.csv("data/expTable.csv",stringsAsFactors = F)

#get the accession numbers which are human
accession.human<-completed[ completed$Species=="Human"&completed$platform!="RNASeq","ID"]
foldChangeOnly.human<-completed[ completed$Species=="Human"&completed$platform!="RNASeq","foldChangeOnly"]
platform.human<-completed[ completed$Species=="Human"&completed$platform!="RNASeq","platform"]

accessionList<-c()
accessionListAll<-c()

getOutputFile<-function(accession,foldChangeOnly){
  print(accession)
  file<-paste0("output/",accession,"-workflowOutput.csv")
  data<-read.csv(file,stringsAsFactors = F)
  downloadUrls<-data$foldChange[-1]
  
  foldChanges <- sapply(downloadUrls,function(x) {
    
    downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",x)
    destFile<-tempfile()
    download.file(downloadUrl,destFile)
    foldChange<-read.delim(destFile)[,1:2]
    colnames(foldChange)[1]<-"GeneSymbol"
    foldChange <- foldChange %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.max(abs(2))) %>% as.data.frame
    foldChange
    
    
  },simplify = F)
  accessionList<<-c(accessionList,rep(accession,length(foldChanges)))
  accessionListAll<<-c(accessionListAll,rep(accession,length(foldChanges)))
  foldChanges <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                        foldChanges)
  
  if (foldChangeOnly==FALSE){
  foldChanges.pval<- sapply(downloadUrls,function(x) {
    
    downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",x)
    destFile<-tempfile()
    download.file(downloadUrl,destFile)
    foldChange<-read.delim(destFile)[,1:3]
    colnames(foldChange)[1]<-"GeneSymbol"
    foldChange <- foldChange %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.min(abs(3))) %>% as.data.frame
    foldChange
    
    
  },simplify = F)
  foldChanges.pval <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                             foldChanges.pval)
  } else{
    foldChanges.pval<-data.frame(GeneSymbol=NA,t(rep(c(NA,NA),length(downloadUrls))))
    
  }
  downloadUrl<-data$chrDirTable[1]
  downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",downloadUrl)
  destFile<-tempfile()
  download.file(downloadUrl,destFile)
  if(file.size(destFile)<5){
    chrDirs <-  vector("list", length(downloadUrls))
  } else{
    chrDir<-read.delim(destFile)
    colnames(chrDir)[1]<-"GeneSymbol"
    chrDirs<- lapply(seq_along(downloadUrls),function(x,chrDir) {
      
      chrDir <- chrDir[ chrDir$comparisonNumber == x,]
      chrDir <- chrDir %>%
        group_by(GeneSymbol) %>%
        dplyr::slice(which.max(abs(2))) %>% as.data.frame
      chrDir
    },chrDir)
  }

  return(list(foldChanges=foldChanges,foldChanges.pval=foldChanges.pval,chrDir=chrDirs))
}

foldChanges<-mapply(FUN = getOutputFile,accession.human,foldChangeOnly.human)
foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.human <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]

foldChangeTable.human<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges)
foldChangeTable.pval.human<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges.pval)
names(chrDirs.human)<- accession.human




#get the accession numbers which are mouse
accession.mouse<-completed[ completed$Species=="Mouse"&completed$platform!="RNASeq","ID"]
foldChangeOnly.mouse<-completed[ completed$Species=="Mouse"&completed$platform!="RNASeq","foldChangeOnly"]

homology<-read.delim("ReferenceData/Homology.mouse.txt")

getOutputFileMouse<-function(accession,foldChangeOnly){
  print(accession)
  file<-list.files(pattern=accession,path = "output/",full.names = T)
  data<-read.csv(file,stringsAsFactors = F)
  downloadUrls<-data$foldChange[-1]
  foldChanges <- sapply(downloadUrls,function(x) {
    
    downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",x)
    destFile<-tempfile()
    download.file(downloadUrl,destFile)
    foldChange<-read.delim(destFile)[,1:2]
    colnames(foldChange)[1]<-"GeneSymbol"
    foldChange<-foldChange[ foldChange$GeneSymbol %in% homology[,3],]
    homoloGene<-homology[match(foldChange$GeneSymbol,homology[,3]),4]
    foldChange$GeneSymbol<-homoloGene
    
    foldChange <- foldChange %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.max(abs(2))) %>% as.data.frame
    foldChange
    
  },simplify = F)
  
  accessionList<<-c(accessionList,rep(accession,length(foldChanges)))
  accessionListAll<<-c(accessionListAll,rep(accession,length(foldChanges)))
  foldChanges <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                        foldChanges)
  
  if (foldChangeOnly==FALSE){
    foldChanges.pval<- sapply(downloadUrls,function(x) {
      
      downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",x)
      destFile<-tempfile()
      download.file(downloadUrl,destFile)
      foldChange<-read.delim(destFile)[,1:3]
      colnames(foldChange)[1]<-"GeneSymbol"
      foldChange<-foldChange[ foldChange$GeneSymbol %in% homology[,3],]
      homoloGene<-homology[match(foldChange$GeneSymbol,homology[,3]),4]
      foldChange$GeneSymbol<-homoloGene
      
      foldChange <- foldChange %>%
        group_by(GeneSymbol) %>%
        dplyr::slice(which.min(abs(3))) %>% as.data.frame
      foldChange
      
    },simplify = F)
    
    foldChanges.pval <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                               foldChanges.pval)
  } else{
    foldChanges.pval<-data.frame(GeneSymbol=NA,t(rep(c(NA,NA),length(downloadUrls))))
    
  }
  downloadUrl<-data$chrDirTable[1]
  downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",downloadUrl)
  destFile<-tempfile()
  download.file(downloadUrl,destFile)
  if(file.size(destFile)<5){
    chrDirs <-  vector("list", length(downloadUrls)) 
  } else{
    chrDir<-read.delim(destFile)
    colnames(chrDir)[1]<-"GeneSymbol"
    chrDirs<- lapply(seq_along(downloadUrls),function(x,chrDir) {

      chrDir <- chrDir[ chrDir$comparisonNumber == x,]
      
      chrDir<-chrDir[ chrDir$GeneSymbol %in% homology[,3],]
      homoloGene<-homology[match(chrDir$GeneSymbol,homology[,3]),4]
      chrDir$GeneSymbol<-homoloGene
      
      chrDir <- chrDir %>%
        group_by(GeneSymbol) %>%
        dplyr::slice(which.max(abs(2))) %>% as.data.frame
      chrDir
    },chrDir)
  }
  return(list(foldChanges=foldChanges,foldChanges.pval=foldChanges.pval,chrDir=chrDirs))
  
}

foldChanges<-mapply(FUN = getOutputFileMouse,accession.mouse,foldChangeOnly.mouse)
foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.mouse <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]

names(chrDirs.mouse)<- accession.mouse

foldChangeTable.mouse<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges)
foldChangeTable.pval.mouse<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges.pval)
foldChangeTable.micro<-merge(foldChangeTable.human,foldChangeTable.mouse, by = "GeneSymbol", all = TRUE)
foldChangeTable.pval.micro<-merge(foldChangeTable.pval.human,foldChangeTable.pval.mouse, by = "GeneSymbol", all = TRUE)

#get the accession numbers which are from rat
accession.rat<-completed[ completed$Species=="Rat"&completed$platform!="RNASeq","ID"]
foldChangeOnly.rat<-as.logical(completed[ completed$Species=="Rat"&completed$platform!="RNASeq","foldChangeOnly"])

homology<-read.delim("ReferenceData/Homology.rat.txt")


foldChanges<-mapply(FUN = getOutputFileMouse,accession.rat,foldChangeOnly.rat)
foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.rat <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]


names(chrDirs.rat)<- accession.rat

foldChangeTable.rat<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                            foldChanges)
foldChangeTable.pval.rat<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                            foldChanges.pval)


foldChangeTable.micro<-merge(foldChangeTable.micro,foldChangeTable.rat,by = "GeneSymbol", all = TRUE)
foldChangeTable.pval.micro<-merge(foldChangeTable.pval.micro,foldChangeTable.pval.rat,by = "GeneSymbol", all = TRUE)


homology<-read.delim("ReferenceData/Homology.pig.txt",stringsAsFactors = F)

#get the accession numbers which are from pig
accession.pig<-completed[ completed$Species=="Pig"&completed$platform!="RNASeq","ID"]
foldChangeOnly.pig<-as.logical(completed[ completed$Species=="Pig"&completed$platform!="RNASeq","foldChangeOnly"])

foldChanges<-mapply(FUN = getOutputFileMouse,accession.pig,foldChangeOnly.pig)
foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.pig <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]

names(chrDirs.pig) <- accession.pig

foldChangeTable.pig <-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                             foldChanges)
foldChangeTable.pval.pig <-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                                  foldChanges.pval)

foldChangeTable.micro<-merge(foldChangeTable.micro,foldChangeTable.pig,by = "GeneSymbol", all = TRUE)
foldChangeTable.pval.micro<-merge(foldChangeTable.pval.micro,foldChangeTable.pval.pig,by = "GeneSymbol", all = TRUE)

##########
#Accession from cow
homology<-read.delim("ReferenceData/Homology.cow.txt",stringsAsFactors = F)

#get the accession numbers which are from pig
accession.cow<-completed[ completed$Species=="Cow"&completed$platform!="RNASeq","ID"]
foldChangeOnly.cow<-as.logical(completed[ completed$Species=="Cow"&completed$platform!="RNASeq","foldChangeOnly"])

foldChanges<-mapply(FUN = getOutputFileMouse,accession.cow,foldChangeOnly.cow)
foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.cow <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]

names(chrDirs.cow) <- accession.cow

foldChangeTable.cow <-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                            foldChanges)
foldChangeTable.pval.cow <-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                                 foldChanges.pval)

foldChangeTable.micro<-merge(foldChangeTable.micro,foldChangeTable.cow,by = "GeneSymbol", all = TRUE)
foldChangeTable.pval.micro<-merge(foldChangeTable.pval.micro,foldChangeTable.pval.cow,by = "GeneSymbol", all = TRUE)



#get the RNA-seq data
accessionList<-c()
#get the accession numbers which are human
accession.human.rna<-completed[ completed$Species=="Human" & completed$platform=="RNASeq","ID"]
foldChangeOnly.human<-completed[ completed$Species=="Human"& completed$platform=="RNASeq","foldChangeOnly"]


getOutputRNAFile<-function(accession,foldChangeOnly){
  print(accession)
  file<-list.files(pattern=paste0(accession,"-workflowOutput.tsv"),path = "output/",full.names = T)
  data<-read.delim(file,stringsAsFactors = F)
  downloadUrls<-data$foldChange[-1]
  foldChanges <- sapply(downloadUrls,function(x) {
    
    downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",x)
    downloadUrl<-gsub(pattern = 'http://localhost:8080/',replacement = "http://localhost:8080/datasets/",downloadUrl)
    destFile<-tempfile()
    download.file(downloadUrl,destFile)
    foldChange<-read.delim(destFile)[,1:2]
    colnames(foldChange)[1]<-"GeneSymbol"
    foldChange <- foldChange %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.max(abs(2))) %>% as.data.frame
    foldChange
    
    
  },simplify = F)
  accessionList<<-c(accessionList,rep(accession,length(foldChanges)))
  accessionListAll<<-c(accessionListAll,rep(accession,length(foldChanges)))
  foldChanges <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                        foldChanges)
  
  if (foldChangeOnly==FALSE){
    foldChanges.pval<- sapply(downloadUrls,function(x) {
      
      downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",x)
      downloadUrl<-gsub(pattern = 'http://localhost:8080/',replacement = "http://localhost:8080/datasets/",downloadUrl)
      destFile<-tempfile()
      download.file(downloadUrl,destFile)
      foldChange<-read.delim(destFile)[,1:3]
      colnames(foldChange)[1]<-"GeneSymbol"
      foldChange <- foldChange %>%
        group_by(GeneSymbol) %>%
        dplyr::slice(which.min(abs(3))) %>% as.data.frame
      foldChange
      
      
    },simplify = F)
    foldChanges.pval <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                               foldChanges.pval)
  } else{
    foldChanges.pval<-data.frame(GeneSymbol=NA,t(rep(c(NA,NA),length(downloadUrls))))
    
  }
  downloadUrl<-data$chrDirTable[1]
  downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",downloadUrl)
  downloadUrl<-gsub(pattern = 'http://localhost:8080/',replacement = "http://localhost:8080/datasets/",downloadUrl)
  destFile<-tempfile()
  download.file(downloadUrl,destFile)
  if(file.size(destFile)<5){
    chrDirs <-  vector("list", length(downloadUrls)) 
  } else{
    chrDir<-read.delim(destFile)
    chrDirs<- lapply(seq_along(downloadUrls),function(x,chrDir) {
    chrDir <- chrDir[ chrDir$comparisonNumber == x,]

            chrDir <- chrDir %>%
        group_by(gene_name) %>%
        dplyr::slice(which.max(abs(2))) %>% as.data.frame
            
      chrDir$ID <- chrDir$gene_name
      chrDir$gene_name = NULL
      chrDir
    },chrDir)
  }
  return(list(foldChanges=foldChanges,foldChanges.pval=foldChanges.pval,chrDir=chrDirs))
}


foldChanges<-mapply(FUN = getOutputRNAFile,accession.human.rna,foldChangeOnly.human)
foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.rnaseq.human <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]

names(chrDirs.rnaseq.human) <- accession.human.rna

foldChangeTable.human<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges)
foldChangeTable.human.pval<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges.pval)



getOutputFileRNAMouse<-function(accession,foldChangeOnly){
  print(accession)
  file<-list.files(pattern=paste0(accession,"-workflowOutput.tsv"),path = "output/",full.names = T)
  data<-read.delim(file,stringsAsFactors = F)
  downloadUrls<-data$foldChange[-1]
  foldChanges <- sapply(downloadUrls,function(x) {
    
    downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",x)
    downloadUrl<-gsub(pattern = 'http://localhost:8080/',replacement = "http://localhost:8080/datasets/",downloadUrl)
    destFile<-tempfile()
    download.file(downloadUrl,destFile)
    foldChange<-read.delim(destFile)[,1:2]
    colnames(foldChange)[1]<-"GeneSymbol"
    foldChange<-foldChange[ foldChange$GeneSymbol %in% homology[,3],]
    homoloGene<-homology[match(foldChange$GeneSymbol,homology[,3]),4]
    foldChange$GeneSymbol<-homoloGene
    
    foldChange <- foldChange %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.max(abs(2))) %>% as.data.frame
    foldChange
    
  },simplify = F)
  
  accessionList<<-c(accessionList,rep(accession,length(foldChanges)))
  accessionListAll<<-c(accessionListAll,rep(accession,length(foldChanges)))
  
  foldChanges <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                        foldChanges)
  
  if (foldChangeOnly==FALSE){
    
  foldChanges.pval <- sapply(downloadUrls,function(x) { 
  
  downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",x)
  downloadUrl<-gsub(pattern = 'http://localhost:8080/',replacement = "http://localhost:8080/datasets/",downloadUrl)
  destFile<-tempfile()
  download.file(downloadUrl,destFile)
  foldChange<-read.delim(destFile)[,1:3]
  colnames(foldChange)[1]<-"GeneSymbol"
  foldChange<-foldChange[ foldChange$GeneSymbol %in% homology[,3],]
  homoloGene<-homology[match(foldChange$GeneSymbol,homology[,3]),4]
  foldChange$GeneSymbol<-homoloGene
  
  foldChange <- foldChange %>%
    group_by(GeneSymbol) %>%
    dplyr::slice(which.min(abs(3))) %>% as.data.frame
  foldChange
  
  },simplify = F)
  
  foldChanges.pval <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                             foldChanges.pval)
} else{
  foldChanges.pval<-data.frame(GeneSymbol=NA,t(rep(c(NA,NA),length(downloadUrls))))
  
}
  downloadUrl<-data[1,"chrDirTable"]
  downloadUrl<-gsub(pattern = '/?preview=True',replacement = "",downloadUrl)
  downloadUrl<-gsub(pattern = 'http://localhost:8080/',replacement = "http://localhost:8080/datasets/",downloadUrl)
  destFile<-tempfile()
  download.file(downloadUrl,destFile)
  if(file.size(destFile)<5){
    chrDirs <-  vector("list", length(downloadUrls)) 
  } else{
    chrDir<-read.delim(destFile)
    chrDirs<- lapply(seq_along(downloadUrls),function(x,chrDir) {
      chrDir <- chrDir[ chrDir$comparisonNumber == x,]
      chrDir<-chrDir[ chrDir$gene_name %in% homology[,3],]
      homoloGene<-homology[match(chrDir$gene_name,homology[,3]),4]
      chrDir$ID<-homoloGene

      chrDir <- chrDir %>%
        group_by(ID) %>%
        dplyr::slice(which.max(abs(2))) %>% as.data.frame
      chrDir$gene_name <- NULL
      chrDir
    },chrDir)
  }
  return(list(foldChanges=foldChanges,foldChanges.pval=foldChanges.pval,chrDir=chrDirs))

}


#get the accession numbers which are mouse
accession.mouse.rna<-completed[ completed$Species=="Mouse"& completed$platform=="RNASeq","ID"]
foldChangeOnly.mouse<-completed[ completed$Species=="Mouse"& completed$platform=="RNASeq","foldChangeOnly"]

homology<-read.delim("ReferenceData/Homology.mouse.txt")

foldChanges<-mapply(FUN = getOutputFileRNAMouse,accession.mouse.rna,foldChangeOnly.mouse)

foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.rnaseq.mouse <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]

names(chrDirs.rnaseq.mouse) <- accession.mouse.rna

foldChangeTable.mouse<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges)

foldChangeTable.pval.mouse<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges.pval)

foldChangeTable.rna<-merge(foldChangeTable.human,foldChangeTable.mouse, by.x= "GeneSymbol",by.y="GeneSymbol",all = TRUE)
foldChangeTable.pval.rna<-merge(foldChangeTable.human.pval,foldChangeTable.pval.mouse, by.x= "GeneSymbol",by.y="GeneSymbol",all = TRUE)

#get the accession numbers which are pig
accession.pig.rna<-completed[ completed$Species=="Pig"& completed$platform=="RNASeq","ID"]
foldChangeOnly.pig<-completed[ completed$Species=="Pig"& completed$platform=="RNASeq","foldChangeOnly"]

homology<-read.delim("ReferenceData/Homology.pig.txt",stringsAsFactors = F)

foldChanges<-mapply(FUN = getOutputFileRNAMouse,accession.pig.rna,foldChangeOnly.pig)

foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.rnaseq.pig <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]

names(chrDirs.rnaseq.pig) <- accession.pig.rna

foldChangeTable.pig<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                              foldChanges)

foldChangeTable.pval.pig<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                                   foldChanges.pval)

foldChangeTable.rna<-merge(foldChangeTable.rna,foldChangeTable.pig, by.x= "GeneSymbol",by.y="GeneSymbol",all = TRUE)
foldChangeTable.pval.rna<-merge(foldChangeTable.pval.rna,foldChangeTable.pval.pig, by.x= "GeneSymbol",by.y="GeneSymbol",all = TRUE)


####
#Rat RNA-seq

#get the accession numbers which are rat
accession.rat.rna<-completed[ completed$Species=="Rat"& completed$platform=="RNASeq","ID"]
foldChangeOnly.rat<-completed[ completed$Species=="Rat"& completed$platform=="RNASeq","foldChangeOnly"]

homology<-read.delim("ReferenceData/Homology.rat.txt",stringsAsFactors = F)

foldChanges<-mapply(FUN = getOutputFileRNAMouse,accession.rat.rna,foldChangeOnly.rat)

foldChanges.pval <- foldChanges[seq(2,length(foldChanges),by=3)]
chrDirs.rnaseq.rat <- foldChanges[seq(3,length(foldChanges),by=3)]
foldChanges <- foldChanges[seq(1,length(foldChanges),by=3)]

names(chrDirs.rnaseq.rat) <- accession.rat.rna

foldChangeTable.rat<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                            foldChanges)

foldChangeTable.pval.rat<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                                 foldChanges.pval)

foldChangeTable.rna<-merge(foldChangeTable.rna,foldChangeTable.rat, by.x= "GeneSymbol",by.y="GeneSymbol",all = TRUE)
foldChangeTable.pval.rna<-merge(foldChangeTable.pval.rna,foldChangeTable.pval.rat, by.x= "GeneSymbol",by.y="GeneSymbol",all = TRUE)

#####################
#merge all the data frames together

foldChangeTable<-merge(foldChangeTable.micro,foldChangeTable.rna, by.y=  "GeneSymbol",by.x="GeneSymbol",all = TRUE)
foldChangeTable.pval<-merge(foldChangeTable.pval.micro,foldChangeTable.pval.rna, by.y=  "GeneSymbol",by.x="GeneSymbol",all = TRUE)
foldChangeTable <- foldChangeTable[ !is.na(foldChangeTable$GeneSymbol),]
foldChangeTable.pval <- foldChangeTable.pval[ !is.na(foldChangeTable.pval$GeneSymbol),]

rownames(foldChangeTable)<-foldChangeTable$GeneSymbol
foldChangeTable$GeneSymbol=NULL

rownames(foldChangeTable.pval)<-foldChangeTable.pval$GeneSymbol
foldChangeTable.pval$GeneSymbol=NULL

xx<-unlist(as.list(org.Hs.egSYMBOL))
foldChangeTable <- foldChangeTable[ rownames(foldChangeTable) %in% xx,]
foldChangeTable.pval <- foldChangeTable.pval[ rownames(foldChangeTable.pval) %in% xx,]

saveRDS(foldChangeTable,file="foldChangeTable.RDS")
saveRDS(foldChangeTable.pval,file="foldChangeTablePVal.RDS")


as.data.frame(accessionListAll) %>%  group_by(accessionListAll) %>%  mutate(myorder = 1:n()) %>% as.data.frame -> accessions
colnames(accessions) <-  c("accession","comparison")
datasetID <- paste(accessions$accession,accessions$comparison,sep="_")



comparisons<- sapply(unique(accessions$accession),function(accession){
  print(accession)
  comparisonsTable <- read.delim(paste0("~/Sybil/Shiny/data/",accession,"/",accession,"_comparisons.txt"))
  sapply(1:nrow(comparisonsTable),function(i) paste0(comparisonsTable[i,"Numerator"],"vs",comparisonsTable[i,"Denominator"]))
})

comparisons <- unlist(comparisons)
accessions$comparisonsText <- comparisons
colnames(foldChangeTable) <- comparisons

write.table(accessions,file="accessions.txt",col.names=T,row.names=F,sep="\t",quote=F)


getZScores<-function(profile,C){
  zscores<-C[-profile,profile]
  rownames(zscores)<-rownames(C[-profile,profile])
  #zscores<-zscores[order(zscores[,1]),]
  return(zscores)
}


cos.sim <- function(ix,X) {
  A = X[ix[1],]
  B = X[ix[2],]
  return( sum(A*B,na.rm = T)/sqrt(sum(A^2,na.rm = T)*sum(B^2,na.rm = T)) )
}   
n <- nrow(t(foldChangeTable)) 
cmb <- expand.grid(i=1:n, j=1:n) 
C <- matrix(apply(cmb,1,cos.sim,t(foldChangeTable)),n,n)
colnames(C)<-colnames(foldChangeTable)
rownames(C)<-colnames(foldChangeTable)

accessions$combined <- paste0(accessions$accession,"_",accessions$comparison)
zscores<-lapply(1:nrow(C),getZScores,C)
names(zscores)<-datasetID
saveRDS(zscores,file="simZScores.RDS")



jaccard<-function(set1,set2){
  I <- length(intersect(set1,set2))
  S <- I/(length(set1)+length(set2)-I)
  return(S)
}

signedJaccard<-function(set1Up,set1Down,set2Up,set2Down){
  signedJaccard<-(jaccard(set1Up,set2Up)+jaccard(set1Down,set2Down)-jaccard(set1Up,set2Down)-jaccard(set1Down,set2Up))/2
  return(signedJaccard)
}

signedJaccardCreate<-function(index,foldChangeTableUp,foldChangeTableDown,foldChangeTable){
  
  i<-index[1]
  j<-index[2]
  
  fc1 <- foldChangeTable[,i,drop=F]
  fc2 <- foldChangeTable[,j,drop=F]
  
  #find the overlap of expressed genes in both datasets
  fc1 <- rownames(fc1[rownames(fc1) %in% rownames(na.omit(fc2)),,drop=F])
  fc2 <-rownames(fc2[rownames(fc2) %in% fc1,,drop=F])
  
  set1Up<-names(foldChangeTableUp[,i][foldChangeTableUp[,i]!=0])
  set2Up<-names(foldChangeTableUp[,j][foldChangeTableUp[,j]!=0])
  set1Down<-names(foldChangeTableDown[,i][foldChangeTableDown[,i]!=0])
  set2Down<-names(foldChangeTableDown[,j][foldChangeTableDown[,j]!=0])
  
  #filter to keep only the genes expressed in both datasets
  set1Up <- set1Up[ set1Up %in% fc1]
  set2Up <- set2Up[ set2Up %in% fc2]
  set1Down <- set1Down[ set1Down %in% fc1]
  set2Down <- set2Down[ set2Down %in% fc2]
  
  #no differentially expressed genes
  if ((length(set1Up) == 0 & length(set1Down) == 0)  | (length(set2Up) == 0 & length(set2Down) == 0) ){
    return(0)
  }
  
  
  return(signedJaccard(set1Up,set1Down,set2Up,set2Down))
}

foldChangeTable <- readRDS("foldChangeTable.RDS")
colnames(foldChangeTable)<-datasetID
foldChangeTableUp<-ifelse(foldChangeTable>=log2(1.5),1,0)
foldChangeTableDown<-ifelse(foldChangeTable<=log2(1/1.5),1,0)

n <- nrow(t(foldChangeTable)) 
cmb <- expand.grid(i=1:n, j=1:n) 
C <- matrix(apply(cmb,1,signedJaccardCreate,foldChangeTableUp,foldChangeTableDown,foldChangeTable),n,n)
colnames(C)<-colnames(foldChangeTable)
rownames(C)<-colnames(foldChangeTable)

zscores<-lapply(1:nrow(C),getZScores,C)
names(zscores)<-datasetID
saveRDS(zscores,file="jaccZScores.RDS")

signedJaccardSigCreate<-function(index,foldChangeListUp,foldChangeListDown,foldChangeList){
  i<-index[1]
  j<-index[2]
  
  #no replicates for differential expression analysis
  if (all(is.na(foldChangeList[[i]])) | all(is.na(foldChangeList[[j]])) ){
    return(NA)
  }
  
  fc1 <- foldChangeList[[i]]
  fc2 <- foldChangeList[[j]]
  
  #find the overlap of expressed genes in both datasets
  fc1 <- fc1[rownames(fc1) %in% rownames(na.omit(fc2)),]
  fc2 <-fc2[rownames(fc2) %in% rownames(na.omit(fc1)),]
  
  #get the up and down genes
  set1Up<-rownames(foldChangeListUp[[i]])
  set2Up<-rownames(foldChangeListUp[[j]])
  set1Down<-rownames(foldChangeListDown[[i]])
  set2Down<-rownames(foldChangeListDown[[j]])
  
  #filter to keep only the genes expressed in both datasets
  set1Up <- set1Up[ set1Up %in% rownames(fc1)]
  set2Up <- set2Up[ set2Up %in% rownames(fc2)]
  set1Down <- set1Down[ set1Down %in% rownames(fc1)]
  set2Down <- set2Down[ set2Down %in% rownames(fc2)]
  
  
  #no differentially expressed genes
  if ((length(set1Up) == 0 & length(set1Down) == 0)  | (length(set2Up) == 0 & length(set2Down) == 0) ){
    return(NA)
  }
  
  
  return(signedJaccard(set1Up,set1Down,set2Up,set2Down))
}

foldChangeList.pval <- lapply(seq(1, ncol(foldChangeTable.pval), by=2), function(i)
foldChangeTable.pval[i: pmin((i+1), ncol(foldChangeTable.pval))])


foldChangeListUp.pval<-lapply(foldChangeList.pval,function(x) na.omit(x[ x[,1] >= log2(1.5) & x[,2] <= 0.05,]))
foldChangeListDown.pval<-lapply(foldChangeList.pval,function(x) na.omit(x[ x[,1] <= log2(1/1.5) & x[,2] <= 0.05,]))

n <- length(foldChangeList.pval)
cmb <- expand.grid(i=1:n, j=1:n)
C <- matrix(apply(cmb,1,signedJaccardSigCreate,foldChangeListUp.pval,foldChangeListDown.pval,foldChangeList.pval),n,n)
colnames(C)<- datasetID
rownames(C)<- datasetID


zscores<-lapply(1:nrow(C),getZScores,C)
names(zscores)<-datasetID

saveRDS(zscores,file="jaccPvalZScores.RDS")

#get all the chdir together and do signed jaccard simimilarity
chrDirs <- c(chrDirs.human,chrDirs.mouse,chrDirs.rat,chrDirs.pig,chrDirs.cow,chrDirs.rnaseq.human,chrDirs.rnaseq.mouse,chrDirs.rnaseq.rat,chrDirs.rnaseq.pig)

chrDirs <- unlist(chrDirs,recursive = F)
names(chrDirs)<-datasetID
save(chrDirs,file="chrDirs.Rdata")
saveRDS(chrDirs,"chrDirs.RDS")



signedJaccardCreateChrDir<-function(index,chrDirs,foldChangeList.pval){
  i<-index[1]
  j<-index[2]
  
  if (is.na( chrDirs[i])| is.na( chrDirs[j])){
    return(NA)
  } 
  
  chrDir1 <- chrDirs[[i]]
  chrDir2 <- chrDirs[[j]]
  
  fc1 <- foldChangeList.pval[[i]]
  fc2 <- foldChangeList.pval[[j]]
  
  #find the overlap of expressed genes in both datasets
  fc1 <- fc1[rownames(fc1) %in% rownames(na.omit(fc2)),]
  fc2 <-fc2[rownames(fc2) %in% rownames(na.omit(fc1)),]
  
  set1Up <- chrDir1[chrDir1$chrDir>0,1]
  set2Up<-chrDir2[chrDir2$chrDir>0,1]
  set1Down<-chrDir1[chrDir1$chrDir<0,1]
  set2Down<-chrDir2[chrDir2$chrDir<0,1]
  
  #filter to keep only the genes expressed in both datasets
  set1Up <- set1Up[ set1Up %in% rownames(fc1)]
  set2Up <- set2Up[ set2Up %in% rownames(fc2)]
  set1Down <- set1Down[ set1Down %in% rownames(fc1)]
  set2Down <- set2Down[ set2Down %in% rownames(fc2)]
  
  
  
  return(signedJaccard(set1Up,set1Down,set2Up,set2Down))
}

#run the chrdir
n <- length(chrDirs)
cmb <- expand.grid(i=1:n, j=1:n)
C <- matrix(apply(cmb,1,signedJaccardCreateChrDir,chrDirs,foldChangeList.pval),n,n)
colnames(C)<- datasetID
rownames(C)<- datasetID
zscores<-lapply(1:nrow(C),getZScores,C)
names(zscores)<-datasetID

saveRDS(C,file="similarityMatrixChrDir.RDS")
saveRDS(zscores,file="chrDirZScores.RDS")


cosineZScores <- readRDS("simZScores.RDS")
jaccardZScores <- readRDS("jaccZScores.RDS")
zscores <- readRDS(file="jaccPvalZScores.RDS")
chrDirZscores <- readRDS("chrDirZScores.RDS")


mergedZscores<-mapply(FUN = function(x,y,z,a) {
  x<-stack(x)
  colnames(x)<-c("cosine","comparisonID")
  y<-stack(y)
  colnames(y)<-c("jaccard","comparisonID")
  z<-stack(z)
  colnames(z)<-c("sigjaccard","comparisonID")
  a<-stack(a)
  colnames(a)<-c("chrDir","comparisonID")
  
  merged<-Reduce(function(dtf1, dtf2) cbind(dtf1, dtf2),
                 list(x,y,z,a))
  merged<-merged[,c(4,3,1,5,7)]
  return(merged)
},cosineZScores,jaccardZScores,zscores,chrDirZscores,SIMPLIFY = F)

save(mergedZscores,file="mergedZscores.RDS")

foldChangeTable$ID <- rownames(foldChangeTable)
write_feather(foldChangeTable, "foldChangeTable.RDS")

#make a folder for the sharedResponses
dataDir<-"~/Sybil/Shiny/data/similarity/"
dir.create(dataDir)

#write each of the dataframes to the folder
mapply(
  write.table,
  x=mergedZscores, file=paste0(dataDir,names(mergedZscores), ".txt"),row.names=FALSE, sep="\t",quote=F)

#create the foldChange.pval lists for the server
foldChangeListUp.pval<-lapply(foldChangeList.pval,function(x) rownames(na.omit(x[ x[,1] >= log2(1.5) & x[,2] <= 0.05,])))
foldChangeListDown.pval<-lapply(foldChangeList.pval,function(x) rownames(na.omit(x[ x[,1] <= log2(1/1.5) & x[,2] <= 0.05,])))
names(foldChangeListUp.pval) <- datasetID
names(foldChangeListDown.pval) <- datasetID

saveRDS(foldChangeListUp.pval,file="foldChangeListUp.pval.RDS",compress = F)
saveRDS(foldChangeListDown.pval,file="foldChangeListDown.pval.RDS",compress = F)


#create a homology mapping file

#mouse
homology.mouse<-read.delim("ReferenceData/Homology.mouse.txt",stringsAsFactors = F)
#rat
homology.rat<-read.delim("ReferenceData/Homology.rat.txt",stringsAsFactors = F)
#cow
homology.cow<-read.delim("ReferenceData/Homology.cow.txt",stringsAsFactors = F)
#pig
homology.pig<-read.delim("ReferenceData/Homology.pig.txt",stringsAsFactors = F)
#add horse
homology.horse<-read.delim("ReferenceData/Homology.horse.txt",stringsAsFactors = F)
#add horse
homology.zebrafish<-read.delim("ReferenceData/Homology.zebrafish.txt",stringsAsFactors = F)

human2otherspecies <-unstack(homology.mouse[,3:4])
human2otherspecies <- lapply(human2otherspecies,function(x) {
  names(x) <- rep("mouse",length(x))
  x
})

rat <-unstack(homology.rat[,3:4])
rat <- lapply(rat,function(x) {
  names(x) <- rep("rat",length(x))
  x
})

human2otherspecies <- merge.list(human2otherspecies,rat)

pig <-unstack(homology.pig[,3:4])
pig <- lapply(pig,function(x) {
  names(x) <- rep("pig",length(x))
  x
})

human2otherspecies <- merge.list(human2otherspecies,pig)

cow <-unstack(homology.cow[,3:4])
cow <- lapply(cow,function(x) {
  names(x) <- rep("cow",length(x))
  x
})

human2otherspecies <- merge.list(human2otherspecies,cow)

horse <-unstack(homology.horse[,3:4])
horse <- lapply(horse,function(x) {
  names(x) <- rep("horse",length(x))
  x
})

human2otherspecies <- merge.list(human2otherspecies,horse)

zebrafish <-unstack(homology.zebrafish[,3:4])
zebrafish <- lapply(zebrafish,function(x) {
  names(x) <- rep("zebrafish",length(x))
  x
})

human2otherspecies <- merge.list(human2otherspecies,zebrafish)


human2otherspecies <- lapply(human2otherspecies,function(x) stack(x))
human2otherspecies <- map_df(human2otherspecies, ~as.data.frame(.x), .id="id")

colnames(human2otherspecies) <- c("humanGene","queryGene","species")


write_feather(human2otherspecies,path ="human2otherspecies.feather")

foldChangeTable.pval <- readRDS("foldChangeTablePVal.RDS")
foldChangeList.pval <- lapply(seq(1, ncol(foldChangeTable.pval), by=2), function(i)
  foldChangeTable.pval[i: pmin((i+1), ncol(foldChangeTable.pval))])
names(foldChangeList.pval)<-accessions$combined

foldChangeList.pval <- lapply(foldChangeList.pval,function(x) x[,2,drop=F])
foldChangeList.pval <- mapply(FUN = function(x,y) {colnames(x) <- y; x},foldChangeList.pval,names(foldChangeList.pval),SIMPLIFY = F)
pvalTable <- bind_cols(foldChangeList.pval)
pvalTable$GeneName <- rownames(foldChangeList.pval[[1]])
write_feather(pvalTable,path = "pvalTable.feather")

