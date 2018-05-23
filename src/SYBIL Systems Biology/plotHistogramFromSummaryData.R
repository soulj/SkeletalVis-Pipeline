#
# Plots histograms
#

suppressPackageStartupMessages(library(RGalaxy))

plotHistogramFromSummaryData <- function(
  meanSummaryData = GalaxyInputFile(required=TRUE),
  sdSummaryData = GalaxyInputFile(required=TRUE),
  outputPlot = GalaxyOutput("histogram", "png")

) {
  meanDataFrame = read.table(meanSummaryData, header = TRUE, sep = '\t', check.names = FALSE)
  groups = as.vector(meanDataFrame[,1])
  measurement = colnames(meanDataFrame)[1]
  meanData = as.matrix(meanDataFrame)[,-1]
  stdData = as.matrix(read.table(sdSummaryData, header = TRUE, sep = '\t', check.names = FALSE))[,-1]
  class(meanData) = "numeric"
  class(stdData) = "numeric"

  png(outputPlot)
  par(mar=c(5.1,4.1,4.1,7.1), xpd=TRUE)
  histogram = barplot(height = meanData, 
                      beside = TRUE, 
                      ylim = c(0, (max(meanData)+max(stdData))+1),
                      col = c('#339933','#ff9933'), 
                      xlab = names(meanData),
                      ylab = measurement)
  histogram = arrows(histogram, meanData-stdData*2, 
                     histogram, meanData+stdData*2, 
                     angle = 90, 
                     length = 0.05, 
                     code = 3)
  legend('topright', legend = groups, fill=c('#339933','#ff9933'), inset=c(-0.2,0))
  dev.off()
}