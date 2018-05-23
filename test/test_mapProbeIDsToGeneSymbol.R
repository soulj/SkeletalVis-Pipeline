options("width"=200)
source('../src/SYBIL Systems Biology/mapProbeIDsToGeneSymbol.R')
  
     inputdata <- 'resources/probe_id_test_data.txt'
   ##   inputdata <- 'resources/MLIIAnalysisSample.txt'
##inputdata <- 'resources/ML2MutantExpressedSample.txt'

t <- 'out_david.txt'
mapProbeIDsToGeneSymbol(
    inputfile1      = inputdata, 
   
    probeColumnIndex= 1, 
    outputfile1     = t)
#message('DAVID query complete.')
x <- read.table(t, header=FALSE, fill=TRUE)
print("PresentGene:")
dim(x)
head(x)
