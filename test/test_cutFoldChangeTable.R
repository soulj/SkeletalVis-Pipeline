source('../src/SYBIL Systems Biology/cutFoldChangeTable.R')


#test fold change only from affy

sourceFile1<-"resources/TestFoldChangeTable.tabular"

t1<-tempfile()

cutFoldChangeTable(foldChangeTable = sourceFile1,
                   comparisonNumber = 1,
                   cutTable = t1)

cutTable <-read.delim(t1)
# Expect dim n, 3 with column names 
# ID
# genotype/variation: Xbp1Cart{delta}Ex2 HZvsgenotype/variation: Wt HZ_logFC 
# genotype/variation: Xbp1Cart{delta}Ex2 HZvsgenotype/variation: Wt HZ_adj.P.Val
dim(cutTable)
head(cutTable)

t1<-tempfile()

cutFoldChangeTable(foldChangeTable = sourceFile1,
                   comparisonNumber = 2,
                   cutTable = t1)
                   
cutTable <-read.delim(t1)
# Expect dim n, 3 with column names 
# ID
# genotype/variation: ColXN617K HZvsgenotype/variation: Wt HZ_logFC 
# genotype/variation: ColXN617K HZvsgenotype/variation: Wt HZ_adj.P.Val
dim(cutTable)
head(cutTable)

sourceFile1<-"resources/TestFoldChangeTableMissingPVal.tabular"

t1<-tempfile()

cutFoldChangeTable(foldChangeTable = sourceFile1,
                   comparisonNumber = 1,
                   cutTable = t1)

cutTable <-read.delim(t1)
# Expect dim n, 2 with column names 
# ID
# genotype/variation: Xbp1Cart{delta}Ex2 HZvsgenotype/variation: Wt HZ
dim(cutTable)
head(cutTable)

t1<-tempfile()

cutFoldChangeTable(foldChangeTable = sourceFile1,
                   comparisonNumber = 2,
                   cutTable = t1)

cutTable <-read.delim(t1)
# Expect dim n, 2 with column names 
# ID
# genotype/variation: ColXN617K HZvsgenotype/variation: Wt HZ
dim(cutTable)
head(cutTable)
