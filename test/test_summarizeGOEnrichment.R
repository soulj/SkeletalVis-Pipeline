
source('../src/SYBIL Systems Biology/summarizeGOEnrichment.R')
# This function will test the function of sumerizeGOEnrichment.
# sent by Galaxy for empty input fields.
t <- tempfile()

summarizeGOEnrichment(
        topicGeneInputfile  = 'resources/GO/topicGenes.txt',
 backgroundSampleInputfile  = 'resources/PCASample/GSM1072318_pes1_KI_OB_I.CEL',
            goTermTypeIndex = 1, 
          goTermOutputfile1 = t
)

w <- read.table(t,header=TRUE,as.is = TRUE, sep="\t")
print('Test received data:')
dim(w)
head(w)
tail(w)