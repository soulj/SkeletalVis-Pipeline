#test download
t<-tempfile()
sourceFile1<-'resources/RNASeq/PRJNA136103_sampleTable.txt'
accessionNumber = "PRJNA136103"

#run script from shell
command<-paste("Rscript","../alt-src/getRNASeqExpressionData/getRNASeqExpressionData.R","--accessionNumber",accessionNumber,"--sampleTable",sourceFile1,"--downloadedFiles",t)

system(command)

file.exists("fastqFiles/SRR096441.fastq.gz")
head(read.delim(t))
