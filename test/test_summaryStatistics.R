source('../src/SYBIL Systems Biology/summaryStatistics.R')

# Note we expect the function to ignore 'none' input, as this is what is
# sent by Galaxy for empty input fields.
t <- "sample-summaryStatistics.tabular"


summaryStatistics(
  inputTable = 'resources/femurLengthTabularSample.tabular',
  summaryStatistic = 'sd',
  outputfile = t
)

print('Summary data created')