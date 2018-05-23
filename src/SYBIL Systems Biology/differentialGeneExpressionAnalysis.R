#
# Differential gene expression analysis using Limma.
#
suppressPackageStartupMessages(library(RGalaxy))

differentialGeneExpressionAnalysis <- function(
  wtfile1=GalaxyInputFile(required=TRUE),
  wtfile2=GalaxyInputFile(required=FALSE),
  wtfile3=GalaxyInputFile(required=FALSE),
  wtfile4=GalaxyInputFile(required=FALSE),
  mufile1=GalaxyInputFile(required=TRUE),
  mufile2=GalaxyInputFile(required=FALSE),
  mufile3=GalaxyInputFile(required=FALSE),
  mufile4=GalaxyInputFile(required=FALSE),
  diffExpressedGenes=GalaxyOutput("differentiallyExpressedGenes", "tabular")
) {
 
  suppressPackageStartupMessages(library(affy))
  suppressPackageStartupMessages(library(gcrma))
  suppressPackageStartupMessages(library(affyQCReport))
  suppressPackageStartupMessages(library(affyPLM))
  suppressPackageStartupMessages(library(limma))
  
  print('Reading affy batch..')
  # get the sizes of the different types of file
  wtInput <- c(wtfile1, wtfile2, wtfile3, wtfile4)
  wtInput <- wtInput[wtInput != 'None']
  muInput <- c(mufile1, mufile2, mufile3, mufile4)
  muInput <- muInput[muInput != 'None']
  cat("Running with ", length(wtInput), " wt files and ", length(muInput), " mu files\n")
  # create the vector for factoring matrix (could well be improved)
  wtBlock <- vector(, length=length(wtInput))
  wtBlock[!wtBlock] <- 1
  muBlock <- vector(, length=length(muInput))
  muBlock[!muBlock] <- 2
  typeVector <- append(wtBlock, muBlock)
  cat("Vector for factoring matrix is ", typeVector, "\n")

  # create the actual input matrix
  affyInput <- c(wtfile1, wtfile2, wtfile3, wtfile4, mufile1, mufile2, mufile3, mufile4)
  affyInput <- affyInput[affyInput != 'None']
  print(affyInput)
  mydata<- ReadAffy(filenames=affyInput)
  
  #normalization
  eset<-rma(mydata)
  featureCount <- dim(eset)[1]
  cat("Feature count is ", featureCount, "\n")
  #2 replicates in each samples
  design <-model.matrix(~0+factor(typeVector));
  colnames(design)<-c("Wt","Mu")
  contrast.matrix<-makeContrasts(MuvsWt=Mu-Wt,levels=design)

  # Now the contrast matrix is combined with the per-probeset linear model fit.
  fit<-lmFit(eset,design)

  #Extract the linear model fit for the contrasts
  fit2<- contrasts.fit(fit, contrast.matrix)
  ebfit2<-eBayes(fit2)

  #result of fold change of 1.5
  probeset.list<-topTable(ebfit2,coef=1,number=featureCount)
  probe1<-rownames(probeset.list)
  write.table(probeset.list, file=diffExpressedGenes, sep = "\t", quote=FALSE, row.names = TRUE, col.names = FALSE)
}