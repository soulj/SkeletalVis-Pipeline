#
# Uses the bioconductor library to perform PCA analysis for given samples.
#
# 
# A expressionset/ affybatch object is expected
# 
# 





suppressPackageStartupMessages(library(RGalaxy))

performPCA <- function(
  inputfile1 = GalaxyInputFile(required=TRUE),
  inputfile2 = GalaxyInputFile(required=FALSE),
  inputfile3 = GalaxyInputFile(required=FALSE),
  inputfile4 = GalaxyInputFile(required=FALSE),
  inputfile5 = GalaxyInputFile(required=FALSE),
  inputfile6 = GalaxyInputFile(required=FALSE),
  outputPlot = GalaxyOutput("plot", "png"),
  inputfile1Name = GalaxyCharacterParam(),
  inputfile2Name = GalaxyCharacterParam(),
  inputfile3Name = GalaxyCharacterParam(),
  inputfile4Name = GalaxyCharacterParam(),
  inputfile5Name = GalaxyCharacterParam(),
  inputfile6Name = GalaxyCharacterParam() 
  ){
  
  suppressPackageStartupMessages(library(affy))
  suppressPackageStartupMessages(library(lattice))
  
  # ReadAffy supports a character vector of files to read into the batch.
  # The input files includeded  four  samples below:
  # 'resources/PCAsample/GSM1072316_pes1_WT_OB_I.CEL'; 
  # 'resources/PCAsample/GSM1072317_pes1_WT_OB_II.CEL';
  # 'resources/PCAsample/GSM1072318_pes1_KI_OB_I.CEL';
  # 'resources/PCAsample/GSM1072319_pes1_KI_OB_II.CEL'.

  print('Reading affy batch..')
  affyInput <- c(inputfile1, inputfile2, inputfile3, inputfile4, inputfile5, inputfile6)
  affyInput <- affyInput[affyInput != 'None']
  print(affyInput)
  nameInput <- c(inputfile1Name, inputfile2Name, inputfile3Name, inputfile4Name, inputfile5Name, inputfile6Name)
  nameInput <- nameInput[1:length(affyInput)]
  dataAffyBatch <- ReadAffy(filenames=affyInput)
  pData(dataAffyBatch)
  print('Performing PCA analysis..')
  
  ## normalize the affy data
  eset<-rma(dataAffyBatch)
  ## generate the sample expression matrix
  exprs<-exprs(eset)
  ## perform the pca
  pca<-prcomp(t(exprs))
  # look at where the variation is across the Principal Components
  summary(pca)
  # #extract the principal components
  pcs<-data.frame(pca$x)
  # add in the annotation data by merging (just for the first 2 principal components)
  pcs<-merge(pcs[,1:2],pData(dataAffyBatch),by=0)
  
  ##plot the PCA by samples
  png(outputPlot)
  out.plot <- xyplot(PC2~PC1, pcs, group=Row.names, 
                     key = simpleKey(nameInput))
  print(out.plot)
  dev.off()
}
