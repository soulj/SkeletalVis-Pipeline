#Test Affy
sourceFile1<-"resources/E-GEOD-8488.rdata"
sourceFile2<-"resources/E-GEOD-8488_comparisons.txt"
platform <- "Affy"
accession <- "E-GEOD-8488"
directory <- tempdir()

command<-paste("Rscript '../alt-src/microarrayQC/microarrayQC.R' --inputfile",sourceFile1,"--comparisonsTable",sourceFile2,"--platform",platform
               ,"--accession",accession,"--directory",directory,sep=" ")
system(command)

#test 1C-Agilent
sourceFile1<-"resources/GSE30628.rdata"
sourceFile2<-"resources/GSE30628_comparisons.txt"
platform <- "1C-Agilent"
accession <- "GSE30628"
directory <- tempdir()

command<-paste("Rscript '../alt-src/microarrayQC/microarrayQC.R' --inputfile",sourceFile1,"--comparisonsTable",sourceFile2,"--platform",platform
               ,"--accession",accession,"--outputDdrectory",directory,sep=" ")
system(command)

#Test Illumina
sourceFile1<-"resources/GSE72261.rdata"
sourceFile2<-"resources/GSE72261_comparisons.txt"
platform <- "Illumina"
accession <- "GSE72261"
directory <- tempdir()

command<-paste("Rscript '../alt-src/microarrayQC/microarrayQC.R' --inputfile",sourceFile1,"--comparisonsTable",sourceFile2,"--platform",platform
               ,"--accession",accession,"--directory",directory,sep=" ")
system(command)
