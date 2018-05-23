source('../src/SYBIL Systems Biology/plotHistogram.R')

# Note we expect the function to ignore 'none' input, as this is what is
# sent by Galaxy for empty input fields.
t <- "sample-plotHistogram.png"


plotHistogram(
  inputfile1 = 'resources/femurLengthSampleData.csv',
  outputPlot = t
)

print('Plot created')