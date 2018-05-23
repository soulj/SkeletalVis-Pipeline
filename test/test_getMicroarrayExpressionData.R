source('../src/SYBIL Systems Biology/getMicroarrayExpressionData.R')

#test agilent
t<-tempfile()
getMicroarrayExpressionData(accessionNumber="GSE30322", platform = "1C-Agilent", expressionData=t, site="GEO",split = F,grepExpression = F)

objectName<-load(t)
print(get(objectName))

#test illumina
t<-tempfile()
getMicroarrayExpressionData(accessionNumber="GSE46750", platform = "Illumina", expressionData=t,site = "GEO", grepExpression = F,split = F,numberLinesSkip = 4)

objectName<-load(t)
print(get(objectName))

#test illumina grep Expression
t<-tempfile()
getMicroarrayExpressionData(accessionNumber="GSE57218", platform = "Illumina", site="GEO",expressionData=t, grepExpression = T,grepString = "expression",split = F)

objectName<-load(t)
print(get(objectName))

#test affy
t<-tempfile()
getMicroarrayExpressionData(accessionNumber="GSE55457", platform = "Affy", expressionData=t, site="GEO",split = F,grepExpression = F,remove = F)

objectName<-load(t)
print(get(objectName))

#test affy from ArrayExpress
t<-tempfile()
getMicroarrayExpressionData(accessionNumber="E-GEOD-8488", platform = "Affy", expressionData=t, site="ArrayExpress",split = T,splitField = "Hybridization.Name",splitSep="_",grepExpression = F,remove=F)

objectName<-load(t)
print(get(objectName))

#test affyST
t<-tempfile()
getMicroarrayExpressionData(accessionNumber="GSE60162", platform = "Affy-ST", expressionData=t, site="GEO",split = F,grepExpression = F,remove=F)

objectName<-load(t)
print(get(objectName))

#test agilent 2C
t<-tempfile()
getMicroarrayExpressionData(accessionNumber="GSE45793", platform = "2C-Agilent", expressionData=t,site="GEO", numberFileRemove = 2,split = F,grepExpression = F)

objectName<-load(t)
print(get(objectName))
