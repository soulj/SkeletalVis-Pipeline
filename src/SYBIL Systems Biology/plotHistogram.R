#
# Plots histograms
#

suppressPackageStartupMessages(library(RGalaxy))

plotHistogram <- function(
  inputfile1 = GalaxyInputFile(required=TRUE),
  outputPlot = GalaxyOutput("plot", "png")

) {
  inputDataFrame = read.csv(inputfile1)
  
  meanData = tapply(inputDataFrame[,4], list(inputDataFrame[,2],inputDataFrame[,1]), mean)
  stdData = tapply(inputDataFrame[,4], list(inputDataFrame[,2],inputDataFrame[,1]), sd)
  
  png(outputPlot)
  par(mar=c(5.1,4.1,4.1,7.1), xpd=TRUE)
  
  histogram = barplot(height = meanData, 
                      beside = TRUE, 
                      ylim = c(0, (max(meanData)+max(stdData))+1),
                      col = c('#339933','#ff9933'), 
                      xlab = names(inputDataFrame)[1],
                      ylab = names(inputDataFrame)[4])
  histogram = arrows(histogram, meanData-stdData*2, 
                     histogram, meanData+stdData*2, 
                     angle = 90, 
                     length = 0.05, 
                     code = 3)
  legend('topright', legend = levels(inputDataFrame[,2]), fill=c('#339933','#ff9933'), inset=c(-0.2,0))
  dev.off()
}