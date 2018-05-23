source('../src/SYBIL Systems Biology/plotHistogramFromSummaryData.R')

# Note we expect the function to ignore 'none' input, as this is what is
# sent by Galaxy for empty input fields.
t <- "sample-plotHistogramFromSummaryData.png"


plotHistogramFromSummaryData(
  meanSummaryData = 'resources/femurLengthMeanSummary.tabular',
  sdSummaryData = 'resources/femurLengthSDSummary.tabular',
  outputPlot = t
)

print('Plot created')