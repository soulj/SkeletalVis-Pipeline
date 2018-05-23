#
# Creates summary statistics
#

suppressPackageStartupMessages(library(RGalaxy))

summaryStatistics <- function(
  inputTable = GalaxyInputFile(required=TRUE),
  summaryStatistic = GalaxySelectParam(c('Mean'='mn', 'Standard Deviation'='sd'), force_select=TRUE, required=TRUE),
  outputfile = GalaxyOutput("summaryStatistics", "tabular")

) {
  inputDataFrame = read.table(inputTable, header = TRUE, sep= '\t', check.names = FALSE)
  
  if (summaryStatistic == 'mn'){
    summaryData = tapply(inputDataFrame[,4], list(inputDataFrame[,2],inputDataFrame[,1]), mean)
  } else if (summaryStatistic == 'sd') {
    summaryData = tapply(inputDataFrame[,4], list(inputDataFrame[,2],inputDataFrame[,1]), sd)
  }
  
  summaryData = cbind(rownames(summaryData), summaryData)
  rownames(summaryData) = NULL
  colnames(summaryData)[1] = names(inputDataFrame)[4]
  
  write.table(summaryData, file=outputfile, sep='\t', quote=FALSE)
}