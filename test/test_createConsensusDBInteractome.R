options("width"=200)
source('../src/SYBIL Systems Biology/createConsensusDBInteractome.R')

sourcefile <- 'resources/ConsensusDBMusSample.txt'

t <- tempfile()
createConsensusDBInteractome(
  consensusDBData = sourcefile,
  interactome     = t)
w <- read.table(t, header=FALSE, fill=TRUE)
print("Result data:")
dim(w)
head(w, n=10)
#tail(w)
#print(w)
if (dim(w)[1] != 8) {
  cat ("ERROR: Expected 8 rows, but got ", dim(w)[1], "\n")
  quit("status = 1")
}
