source('../src/SYBIL Systems Biology/characteristicDirectionFromArrayExpress.R')


#test fold change only from affy

sourceFile1<-"resources/E-GEOD-8488.rdata"
sourceFile2<-"resources/E-GEOD-8488_comparisons.txt"

t1<-tempfile()

characteristicDirectionFromArrayExpress(inputfile = sourceFile1,
                                       comparisonsTable = sourceFile2,
                                       platform = "Affy",
                                       annotationFile="mouse4302.db",
                                       foldChangeOnly = T,
                                       chrDirTable = t1)



#test p-value calculation and multiple conditions from agilent
sourceFile1<-"resources/GSE30628.rdata"
sourceFile2<-"resources/GSE30628_comparisons.txt"

t1<-tempfile()

characteristicDirectionFromArrayExpress(inputfile = sourceFile1,
                                       comparisonsTable = sourceFile2,
                                       platform = "1C-Agilent",
                                       annotationFile = "mgug4122a.db",
                                       foldChangeOnly = F,
                                       chrDirTable = t1)

chrDirTable <-read.delim(t1)
dim(chrDirTable)
head(chrDirTable)

#test p-value calculation and multiple conditions from Illumina
sourceFile1<-"resources/GSE72261.rdata"
sourceFile2<-"resources/GSE72261_comparisons.txt"

t1<-tempfile()

characteristicDirectionFromArrayExpress(inputfile = sourceFile1,
                                       comparisonsTable = sourceFile2,
                                       platform = "Illumina",
                                       annotationFile="lumiMouseIDMapping",
                                       foldChangeOnly = F,
                                       offset = 0,
                                       chrDirTable = t1)

chrDirTable <-read.delim(t1)
dim(chrDirTable)
head(chrDirTable)


#test p-value calculation and multiple factors that need to be combined (tissue and genotype).
sourceFile1<-"resources/E-GEOD-18647.rdata"
sourceFile2<-"resources/E-GEOD-18647_comparisons.txt"

t1<-tempfile()

characteristicDirectionFromArrayExpress(inputfile = sourceFile1,
                                       comparisonsTable = sourceFile2,
                                       platform = "Affy",
                                       annotationFile="mouse4302.db",
                                       foldChangeOnly = F,
                                       chrDirTable = t1)

chrDirTable <-read.delim(t1)
dim(chrDirTable)
head(chrDirTable)

