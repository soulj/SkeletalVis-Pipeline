
source('../src/SYBIL Systems Biology/plotVennDiagram.R')

# Note we expect the function to compare at least two data sets.
# Note: the file 'resources/ModuleNodes/empty.txt' is for the situation that the algorithms could not detect any nodes in the module.
# We created the empty file which only contains a null value in the table.
t <- "sample-plotVennDigram.pdf"
modules <- tempfile()

plotVennDiagram (
  inputfile1 = 'resources/ModuleNodes/nodeOfBioNetOC.txt',
  inputfile2 = 'resources/ModuleNodes/nodeOfClustExOC.txt',
  inputfile3 = 'resources/ModuleNodes/empty.txt',
  inputfile4 = 'resources/ModuleNodes/empty.txt',
  inputfile5 = 'resources/ModuleNodes/empty.txt',
  intersectionThreshold = 2,
  outputfile1      = modules,
  outputPlot = t
  
)
nodeTable <- read.table(modules, header=FALSE, fill=TRUE)
print("Nodes:")
dim(nodeTable)
head(nodeTable)
print('Plot created')
