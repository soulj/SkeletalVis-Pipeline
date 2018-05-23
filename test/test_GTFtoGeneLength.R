options("width"=200)
source('../src/SYBIL Systems Biology/GTFtoGeneLength.R')
sourcefile1 <- 'resources/RNASeq/humanGTFSample.gtf'



geneExonLengths <- tempfile()

GTFtoGeneLength(gtffile = sourcefile1,species = "Human",geneExonLengths = geneExonLengths)


sessionInfo()

geneExonLengths  <- read.delim(geneExonLengths)
print("geneExonLengths : ")
dim(geneExonLengths)
head(geneExonLengths)
