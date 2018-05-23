
source('../src/SYBIL Systems Biology/extractAllExpressedProbes.R')

# Note we expect the function to ignore 'none' input, as this is what is
# sent by Galaxy for empty input fields.
t <- tempfile()
extractAllExpressedProbes(
  
  ##inputfile1  = 'resources/GPL1261_ch1/GSM469874.CEL.gz',
  ##inputfile2  = 'resources/GPL1261_ch1/GSM469875.CEL.gz',
  inputfile1  = 'resources/GPL1261_ch1/GSM913972_WT1.CEL.gz',
  inputfile2  = 'resources/GPL1261_ch1/GSM913973_WT2.CEL.gz',
  ##inputfile3  = 'None',
  outputfile1 = t
)

w <- read.table(t)
print('Test received data:')
dim(w)
head(w)
tail(w)

t <- tempfile()
extractAllExpressedProbes(
  inputfile1  = 'resources/GPL1261_ch1/GSM549992.CEL.gz',
  ##inputfile1  = 'resources/GPL1261_ch1/GSM469874.CEL.gz',
  ##inputfile2  = 'resources/GPL1261_ch1/GSM469875.CEL.gz',
  ##inputfile1  = 'resources/GPL1261_ch1/GSM913972_WT1.CEL.gz',
  ##inputfile2  = 'resources/GPL1261_ch1/GSM913973_WT2.CEL.gz',
  ##inputfile3  = 'None',
  outputfile1 = t
)

w <- read.table(t)
print('Test received data:')
dim(w)
head(w)
tail(w)