#
# Uses the bioconductor library to read array files and report on protein interactions.
#
# Multiple files in the AffyMetrix batch are supported, only the first is required.
#
# See https://www.bioconductor.org/packages/release/bioc/manuals/affy/man/affy.pdf
# CEL format: http://www.affymetrix.com/estore/support/developer/powertools/changelog/gcos-agcc/cel.html.affx
#

suppressPackageStartupMessages(library(RGalaxy))

extractAllExpressedProbes <- function(
  inputfile1=GalaxyInputFile(required=TRUE),
  inputfile2=GalaxyInputFile(required=FALSE),
  inputfile3=GalaxyInputFile(required=FALSE),
  inputfile4=GalaxyInputFile(required=FALSE),
  inputfile5=GalaxyInputFile(required=FALSE),
  inputfile6=GalaxyInputFile(required=FALSE),
  outputfile1=GalaxyOutput("output", "tabular")
) {
  suppressPackageStartupMessages(library(affy))

  # ReadAffy supports a character vector of files to read into the batch.
  # Note we must remove any input which hasn't been set. This requires some arbitrary
  # knowledge of what Galaxy sends.
  ## The input file 1 includeded  two  samples below:
  ## 'resources/GPL1261_ch1/GSM469874.CEL.gz' and 'resources/GPL1261_ch1/GSM469875.CEL.gz', which were already stored in Galaxy Dataset
 
  print('Reading affy batch..')
  affyInput <- c(inputfile1, inputfile2, inputfile3, inputfile4, inputfile5, inputfile6)
  affyInput <- affyInput[affyInput != 'None']
  print(affyInput)
  dataAffyBatch <- ReadAffy(filenames=affyInput)

  print('Performing MAS 5.0 absolute detection..')
  mas5(dataAffyBatch)
  dataMas5 <-mas5calls(dataAffyBatch)

  # Note that Galaxy can't handle the column names.
  # row names ARE required since this contains the expression name.
  print('Writing output..')
  ## get the table of mas5calls result as follow:
  
  ##                     1415670_at  P  P
  ##   AFFX-TransRecMur/X57349_M_at  P  A
  ##                AFFX-TrpnX-5_at  A  A
  ##  
  
  ## The row name of the table shows the probe id. 
  ## Column 1 and column 2 of the input file 2  indicats whether the transcript was present (P), absent (A), or marginal (M). 
  data1<-assayDataElement(dataMas5,"exprs")
  ##check if any element in a row  of data1 contains "P". If it contained "P", the probe id is present.
  t0<-apply(data1, 1, function(r) any(r %in% c("P")))
  ## remove all probes whose row does not contains "P".
  t1<-as.matrix(t0[t0=="TRUE"])
  presentProbe<-rownames(t1)
  probeid<-as.character(presentProbe)
  presentProbe1<-unique(probeid)
 
  
  ## The output of above eample should be a list of probe ids as follow:
  
  ##                     1415670_at	
  ##   AFFX-TransRecMur/X57349_M_at	
  write.table(probeid, file=outputfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, na="", sep="\t") 
  print(presentProbe1)
  
  
}
