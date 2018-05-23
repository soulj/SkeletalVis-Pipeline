options("width"=200)
source('../src/SYBIL Systems Biology/differentialGeneExpressionAnalysis.R')

t <- tempfile()
differentialGeneExpressionAnalysis(
  wtfile1 = 'resources/Cant1KOMice/1_WT68_(Mouse430_2).CEL',
  wtfile2 = 'resources/Cant1KOMice/2_WT70_(Mouse430_2).CEL',
  wtfile3 = 'resources/Cant1KOMice/3_WT149_(Mouse430_2).CEL',
  wtfile4 = 'None',
  mufile1 = 'resources/Cant1KOMice/4_KO59_(Mouse430_2).CEL',
  mufile2 = 'resources/Cant1KOMice/5_KO62_(Mouse430_2).CEL',
  mufile3 = 'resources/Cant1KOMice/6_KO150_(Mouse430_2).CEL',
  mufile4 = 'None',
  t)
w <- read.table(t, header=FALSE, fill=TRUE)
print("Result data:")
dim(w)
head(w)

t <- tempfile()
differentialGeneExpressionAnalysis(
  wtfile1 = '../../sample-data/ML2\ disease/GSM1072316_pes1_WT_OB_I.CEL',
  wtfile2 = '../../sample-data/ML2\ disease/GSM1072317_pes1_WT_OB_II.CEL',
  wtfile3 = 'None',
  wtfile4 = 'None',
  mufile1 = '../../sample-data/ML2\ disease/GSM1072318_pes1_KI_OB_I.CEL',
  mufile2 = '../../sample-data/ML2\ disease/GSM1072319_pes1_KI_OB_II.CEL',
  mufile3 = 'None',
  mufile4 = 'None',
  t)
w <- read.table(t, header=FALSE, fill=TRUE)
print("Result data:")
dim(w)
head(w)
