
source('../src/SYBIL Systems Biology/performPCA.R')

# Note we expect the function to ignore 'none' input, as this is what is
# sent by Galaxy for empty input fields.
t <- "sample-performPCA.png"

# Note that meta data for file names is passed in as a separate argument.
# This is work in progress to resolve issue 10789.
performPCA(
  inputfile1 = 'resources/PCAsample/GSM1072316_pes1_WT_OB_I.CEL',
  inputfile2 = 'resources/PCAsample/GSM1072317_pes1_WT_OB_II.CEL',
  inputfile3 = 'resources/PCAsample/GSM1072318_pes1_KI_OB_I.CEL',
  inputfile4 = 'resources/PCAsample/GSM1072319_pes1_KI_OB_II.CEL',
  inputfile5 = 'None',
  inputfile6 = 'None',
  outputPlot = t,
  inputfile1Name = "file1Name",
  inputfile2Name = "file2Name",
  inputfile3Name = "file3Name",
  inputfile4Name = "file4Name"
)

print('Plot created')
