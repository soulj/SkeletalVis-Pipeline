# Takes a table of fold changes for multiple comparisons and will return table of fold changes for one comparison

suppressPackageStartupMessages(library(RGalaxy))

cutFoldChangeTable <- function(
  foldChangeTable = GalaxyInputFile(required=TRUE, formatFilter = "tabular"), 
  comparisonNumber = GalaxyNumericParam(required=TRUE),
  cutTable = GalaxyOutput("cutTable","tabular")) {
  
  # Get specified number of columns from table starting at specified index and the first column.
  # 
  # param FCTable The fold change table to get columns from.
  # param columnIndex The start index of the columns you want to get
  # param numberOfColumns Number of columns to return
  getColumns <- function(FCTable, columnIndex, numberOfColumns) {
    cbind(FCTable[1], FCTable[seq(columnIndex, columnIndex+numberOfColumns-1)])
  }
  
  foldChangeTable = read.table(foldChangeTable, header = TRUE, sep = '\t')
  columnNames = names(foldChangeTable)
  if (length(grep('P.Val|padj', columnNames)) == 0) {
    # Only have fold change
    returnTable = getColumns(foldChangeTable, comparisonNumber+1, 1)
  } else {
    returnTable = getColumns(foldChangeTable, comparisonNumber*2, 2)
  }
  write.table(returnTable, file = cutTable, col.names=T, row.names=F, sep="\t", quote=F)
}