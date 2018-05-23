options("width"=200)
source('../src/SYBIL Systems Biology/createSTRINGInteractome.R')

t <- tempfile()
createSTRINGInteractome(
  stringData           = 'resources/10090.protein.actions.test.txt',
 # stringProteinAliases = 'resources/StringMusFullGeneNameList0.txt',
  stringProteinAliases = 'resources/STRINGv10MusFullGeneNameList0.txt',
  confidenceThreashold = 0.7,
  interactome          = t)
w <- read.table(t, header=FALSE, fill=TRUE)
print("Result data:")
dim(w)
head(w)
tail(w)
