source('../src/SYBIL Systems Biology/differentialExpressionFromArrayExpress.R')


#test fold change only from affy

sourceFile1<-"resources/E-GEOD-8488.rdata"
sourceFile2<-"resources/E-GEOD-8488_comparisons.txt"

t1<-tempfile()

differentialExpressionFromArrayExpress(inputfile = sourceFile1,
                                       comparisonsTable = sourceFile2,
                                       platform = "Affy",
                                       annotationFile="mouse4302.db",
                                       foldChangeOnly = T,
                                       foldChangeTable = t1)

foldChangeTable <-read.delim(t1)
dim(foldChangeTable)
head(foldChangeTable)


#test p-value calculation and multiple conditions from agilent
sourceFile1<-"resources/GSE30628.rdata"
sourceFile2<-"resources/GSE30628_comparisons.txt"

t1<-tempfile()

differentialExpressionFromArrayExpress(inputfile = sourceFile1,
                                       comparisonsTable = sourceFile2,
                                       platform = "1C-Agilent",
                                       annotationFile = "mgug4122a.db",
                                       foldChangeOnly = F,
                                       foldChangeTable = t1)

foldChangeTable <-read.delim(t1)
dim(foldChangeTable)
head(foldChangeTable)

#test p-value calculation and multiple conditions from Illumina
sourceFile1<-"resources/GSE72261.rdata"
sourceFile2<-"resources/GSE72261_comparisons.txt"

t1<-tempfile()

differentialExpressionFromArrayExpress(inputfile = sourceFile1,
                                       comparisonsTable = sourceFile2,
                                       platform = "Illumina",
                                       annotationFile="lumiMouseIDMapping",
                                       foldChangeOnly = F,
                                       offset = 0,
                                       foldChangeTable = t1)

foldChangeTable <-read.delim(t1)
dim(foldChangeTable)
head(foldChangeTable)


#test p-value calculation and multiple factors that need to be combined (tissue and genotype).
sourceFile1<-"resources/E-GEOD-18647.rdata"
sourceFile2<-"resources/E-GEOD-18647_comparisons.txt"

t1<-tempfile()

differentialExpressionFromArrayExpress(inputfile = sourceFile1,
                                       comparisonsTable = sourceFile2,
                                       platform = "Affy",
                                       annotationFile="mouse4302.db",
                                       foldChangeOnly = F,
                                       foldChangeTable = t1)

foldChangeTable <-read.delim(t1)
dim(foldChangeTable)
head(foldChangeTable)

