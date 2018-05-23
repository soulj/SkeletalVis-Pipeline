#
# Uses the bioconductor VennDigram library to plot Venn Diagram to show comparison results between given lists of input nodes.
#
# Multiple text files in the text formt were read as input
## This function is to compare genes detected by different functional module detection algorithms.
## It will export the intersection module over certain predefined number of algorithms and show a Venn diagram for the common area distribution.
## We test the function intersection module detected by over 2 different algorithms.
## It is assumed that at least two valid data sets should be imported for comparison.
## 
## Note:
## This script only works if the all the inputs are set. If there is no data for the required slot, then use an 
## empty text file
##

suppressPackageStartupMessages(library(RGalaxy))

plotVennDiagram <- function(
  inputfile1 = GalaxyInputFile(required=TRUE),
  inputfile2 = GalaxyInputFile(required=TRUE),
  inputfile3 = GalaxyInputFile(required=TRUE),
  inputfile4 = GalaxyInputFile(required=TRUE),
  inputfile5 = GalaxyInputFile(required=TRUE),
  intersectionThreshold = GalaxyNumericParam(1, required=TRUE, min=1),
  outputfile1           = GalaxyOutput("IntersectionModuleNode", "tabular"),
  outputPlot = GalaxyOutput("plot", "pdf")
  
) {
  suppressPackageStartupMessages(library(affy))
  suppressPackageStartupMessages(library(lattice))
  suppressPackageStartupMessages(library(VennDiagram))
 
  print('Reading input data..')
  
  ## read the input data set for comparison
  ## Each inputfile is a list of genes which are present in the functional modules that detected by different algorithms.
  ## We want to compare genes detected by different alogtirhms and find the common genes.
  
  input1 <- read.table(inputfile1, sep='\t', as.is=TRUE)
  ## The input of data is a list of gene symbols as follow:
  
  ##  RBP1
  ##  IRF7
  ##  AKT1
  ##  SMAD6
  ##  UBC
  ##  ATM
  
  ## The format of input2 is the same as the input1 data.
  input2 <- read.table(inputfile2, sep='\t', as.is=TRUE)
  ## The format of input3 is the same as the input1 data.
  input3 <- read.table(inputfile3, sep='\t', as.is=TRUE)
  ## The format of input4 is the same as the input1 data.
  input4 <- read.table(inputfile4, sep='\t', as.is=TRUE)
  ## The format of input5 is the same as the input1 data.
  input5 <- read.table(inputfile5, sep='\t', as.is=TRUE)
  
   ## remove duplicate genes in the list
  data1<-as.character(unique(input1[,1]))
  
  ## remove duplicate genes in the list before comparison.
  data2<-as.character(unique(input2[,1]))
  
  ## remove duplicate genes in the list before comparison.
  data3<-as.character(unique(input3[,1]))
 
  ## remove duplicate genes in the list before comparison.
  data4<-as.character(unique(input4[,1]))
  
  ## remove duplicate genes in the list before comparison.
  data5<-as.character(unique(input5[,1]))
  
  ## Start to plot the Venn diagram according to the input data.
  ## Various functions were called to plot the Venn diagram according to the number of valid input data.
  
  print('plot Venn Digram..')
   if(length(data5)>1){

   ## All these five inputfiles are not empty. 
   ## We plot a Venn digaram to compare these five data sets (inputfile1, inputfile2,inputfile3, inputfile4 and inputfile5).
   ##  plot the Venn diagram
   pdf(outputPlot)
   out.plot <-draw.quintuple.venn(
     area1 = length(data1),
     area2 = length(data2), 
     area3 = length(data3),
     area4 = length(data4),
     area5 = length(data5),
     n12=length(intersect(data1,data2)), 
     n13=length(intersect(data1,data3)),
     n14=length(intersect(data1,data4)),
     n15=length(intersect(data1,data5)),
     n23=length(intersect(data2,data3)),
     n24=length(intersect(data2,data4)),
     n25=length(intersect(data2,data5)),
     n34=length(intersect(data3,data4)),
     n35=length(intersect(data3,data5)),
     n45=length(intersect(data4,data5)),
     n123=length(Reduce(intersect, list(data1,data2,data3))),
     n124=length(Reduce(intersect, list(data1,data2,data4))),
     n125=length(Reduce(intersect, list(data1,data2,data5))),
     n134=length(Reduce(intersect, list(data1,data3,data4))),
     n135=length(Reduce(intersect, list(data1,data3,data5))),
     n145=length(Reduce(intersect, list(data1,data4,data5))),
     n234=length(Reduce(intersect, list(data2,data3,data4))),
     n235=length(Reduce(intersect, list(data2,data3,data5))),
     n245=length(Reduce(intersect, list(data2,data4,data5))),
     n345=length(Reduce(intersect, list(data3,data4,data5))),
     n1234=length(Reduce(intersect, list(data1,data2,data3,data4))),
     n1235=length(Reduce(intersect, list(data1,data2,data3,data5))),
     n1245=length(Reduce(intersect, list(data1,data2,data4,data5))),
     n1345=length(Reduce(intersect, list(data1,data3,data4,data5))),
     n2345=length(Reduce(intersect, list(data2,data3,data4,data5))),
     n12345=length(Reduce(intersect, list(data1,data2,data3,data4,data5))),
     category = c("input1", "input2","input3","input4","input5"),
     lwd=rep(2,5),
     lty=rep("solid",5),
     fill = c("dodgerblue","goldenrod2","darkorange2", "red3","orchid3"),
     cat.col =c("dodgerblue","goldenrod2","darkorange2", "red3","orchid3"),
     cat.cex=1.2,
     margin=0.05,
     cex =1,
     ind=TRUE)
   print(out.plot)
   dev.off()
   
   ## export the intersection module detected over certain alogithms, the number of which is defined by the intersection threshold.
   
   ## get common area of certain input data set in the Venn Digram.
   ## common area between genes from the inputfile1 and inputfile2.
   node12<-as.character(intersect(data1,data2))
   ## common area between genes from the inputfile1 and inputfile3.
   node13<-as.character(intersect(data1,data3))
   ## common area between genes from the inputfile1 and inputfile4.
   node14<-as.character(intersect(data1,data4))
   ## common area between genes from the inputfile1 and inputfile5.
   node15<-as.character(intersect(data1,data5))
   ## common area between genes from the inputfile2 and inputfile3.
   node23<-as.character(intersect(data2,data3))
   ## common area between genes from the inputfile2 and inputfile4.
   node24<-as.character(intersect(data2,data4))
   ## common area between genes from the inputfile2 and inputfile5.
   node25<-as.character(intersect(data2,data5))
   ## common area between genes from the inputfile3 and inputfile4.
   node34<-as.character(intersect(data3,data4))
   ## common area between genes from the inputfile3 and inputfile5.
   node35<-as.character(intersect(data3,data5))
   ## common area between genes from the inputfile4 and inputfile5.
   node45<-as.character(intersect(data4,data5))
   ## common area among genes from the inputfile1, inputfile2 and inputfile3.
   node123<-as.character(Reduce(intersect, list(data1,data2,data3)))
   ## common area among genes from the inputfile1, inputfile2 and inputfile4.
   node124<-as.character(Reduce(intersect, list(data1,data2,data4)))
   ## common area among genes from the inputfile1, inputfile2 and inputfile5.
   node125<-as.character(Reduce(intersect, list(data1,data2,data5)))
   ## common area among genes from the inputfile1, inputfile3 and inputfile4.
   node134<-as.character(Reduce(intersect, list(data1,data3,data4)))
   ## common area among genes from the inputfile1, inputfile3 and inputfile5.
   node135<-as.character(Reduce(intersect, list(data1,data3,data4)))
   ## common area among genes from the inputfile1, inputfile4 and inputfile5.
   node145<-as.character(Reduce(intersect, list(data1,data4,data5)))
   ## common area among genes from the inputfile2, inputfile3 and inputfile4.
   node234<-as.character(Reduce(intersect, list(data2,data3,data4)))
   ## common area among genes from the inputfile2, inputfile4 and inputfile5.
   node245<-as.character(Reduce(intersect, list(data2,data4,data5)))
   ## common area among genes from the inputfile2, inputfile3 and inputfile5.
   node235<-as.character(Reduce(intersect, list(data2,data3,data5)))
   ## common area among genes from the inputfile3, inputfile4 and inputfile5.
   node345<-as.character(Reduce(intersect, list(data3,data4,data5)))
   ## common area among genes from the inputfile1, inputfile2, inputfile3 and inputfile4.
   node1234<-as.character(Reduce(intersect, list(data1,data2,data3,data4)))
   ## common area among genes from the inputfile1, inputfile2, inputfile3 and inputfile5.
   node1235<-as.character(Reduce(intersect, list(data1,data2,data3,data5)))
   ## common area among genes from the inputfile1, inputfile2, inputfile4 and inputfile4.
   node1245<-as.character(Reduce(intersect, list(data1,data2,data4,data5)))
   ## common area among genes from the inputfile1, inputfile3, inputfile4 and inputfile4.
   node1345<-as.character(Reduce(intersect, list(data1,data3,data4,data5)))
   ## common area among genes from the inputfile2, inputfile3, inputfile4 and inputfile5.
   node2345<-as.character(Reduce(intersect, list(data2,data3,data4,data5)))
   ## common area among genes from the inputfile1, inputfile2, inputfile3  inputfile4 and inputfile5.
   node12345<-as.character(Reduce(intersect, list(data1,data2,data3,data4,data5)))
   
   ## detect the intersection module detected by at least one alogrithm, which is a combination of all nodes detected by these five algorithms.
   module_over_one=unique(c(data1,data2,data3,data4,data5))
  
   ## detect the intersection module detected by over any two alogrithms, which contains the intersection module detected by any two algortihsm, any three algorithms, any four algorithms or  five algorithms.
   module_over_two=unique(c(node12,node13,node14,node15,node23,node24,node25,node34,node35,node45,node123,node124,node125,node134,node135,node145,node234,node245,node235,node345,node1234,node1235,node1245,node1345,node2345,node12345))
   ## detect the intersection module detected by over any three algorithms, which contains the intersection module detected by any three algorithms,  four algorithms or  five algorithms.
   module_over_three=unique(c(node123,node124,node125,node134,node135,node145,node234,node245,node235,node345,node1234,node1235,node1245,node1345,node2345,node12345))
   
   ## detect the intersection module detected by over any four algorithms, which includes the intersection module is common in any four datasets or all five algorithms.
   module_over_four=unique(c(node1234,node1235,node1245,node1345,node2345,node12345))
   ## detect the intersection module detected by all five datasets
   module_over_five=unique(node12345)
   ## get the intersection module detected by over certain number of algorithms,  which is defined by the intersectionThreshold.
   
   ## In the test example, we want to find the module detected over any four algorithms
   intersection_module5=switch(intersectionThreshold,module_over_one,module_over_two,module_over_three,module_over_four,module_over_five )
   
   ##  The output of the intersection module is a list of genes  which are present in the intersection module as follow:
   
   ##  RBP1
   ##  IRF7
   ##  ID4
   ##  PTGIS
   ##  COL2A1
   ##  PTX3
   write.table(intersection_module5, file=outputfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
 }
  else if(length(data4)>1){
 
    ## The inputfile5 data is empty, but the data in the other four inputfiles are not.
    ## So we plot a Venn digaram to compare four data sets.
    
    ## plot the Venn Diagram
    pdf(outputPlot)
out.plot<-draw.quad.venn(
      area1 = length(data1),
      area2 = length(data2), 
      area3 = length(data3),
      area4 = length(data4),
      n12=length(intersect(data1,data2)), 
      n13=length(intersect(data1,data3)),
      n14=length(intersect(data1,data4)),
      n23=length(intersect(data2,data3)),
      n24=length(intersect(data2,data4)),
      n34=length(intersect(data3,data4)),
      n123=length(Reduce(intersect, list(data1,data2,data3))),
      n124=length(Reduce(intersect, list(data1,data2,data4))),
      n134=length(Reduce(intersect, list(data1,data3,data4))),
      n234=length(Reduce(intersect, list(data2,data3,data4))),
      n1234=length(Reduce(intersect, list(data1,data2,data3,data4))),
      category = c("input1", "input2","input3","input4"),
      fill = c("orange","red", "green","blue"),
      lty="solid",
      cex =2, 
      cat.cet=2,
      cat.col =c("orange","red", "green","blue")
    )
   print(out.plot)
   dev.off()
  ## export the intersection module detected over the intersection threshold.
  ## get common area of certain input data set in the Venn Digram.
  ## common area between genes from the inputfile1 and inputfile2. 
  node12<-as.character(intersect(data1,data2))
  ## common area between genes from the inputfile1 and inputfile3.
  node13<-as.character(intersect(data1,data3))
  ## common area between genes from the inputfile1 and inputfile4.
  node14<-as.character(intersect(data1,data4))
  ## common area between genes from the inputfile2 and inputfile3.
  node23<-as.character(intersect(data2,data3))
  ## common area between genes from the inputfile2 and inputfile4.
  node24<-as.character(intersect(data2,data4))
  ## common area between genes from the inputfile3 and inputfile4.
  node34<-as.character(intersect(data3,data4))
  ## common area among genes from the inputfile1, inputfile2 and inputfile3.
  node123<-as.character(Reduce(intersect, list(data1,data2,data3)))
  ## common area among genes from the inputfile1, inputfile2 and inputfile4.
  node124<-as.character(Reduce(intersect, list(data1,data2,data4)))
  ## common area among genes from the inputfile1, inputfile3 and inputfile4.
  node134<-as.character(Reduce(intersect, list(data1,data3,data4)))
  ## common area among genes from the inputfile2, inputfile3 and inputfile4.
  node234<-as.character(Reduce(intersect, list(data2,data3,data4)))
  ## common area among genes from the inputfile1, inputfile2, inputfile3 and inputfile4. 
  node1234<-as.character(Reduce(intersect, list(data1,data2,data3,data4)))
  
  ## detect the intersection module detected by all four algorithms.
  module_over_one=unique(c(data1,data2,data3,data4))
  ## detect the intersection module detected by over any two algorithms, which contains the intersection module detected by any two algorithms, any three algorithms, or any four algorithms.
  module_over_two<-unique(c(node12, node13, node14, node23, node24, node34, node123, node124,node134, node234, node1234))
  ## detect the intersection module detected by over any three algorithms, which contains the intersection module detected by any three algorithms, or any four algorithms.
  module_over_three<-unique(c(node123,node124,node134,node234,node1234))
  ## detect the intersection module detected by all four algorithms.
  module_over_four<-unique(c(node1234))
  ## get the intersection module detected by over certain number of algorithms,  which is defined by the intersectionThreshold.
  intersection_module4=switch(intersectionThreshold,module_over_one,module_over_two,module_over_three,module_over_four)
  ##  The output of the intersection module is a list of genes  which are present in the intersection module as follow:
  
  ##  RBP1
  ##  IRF7
  ##  ID4
  ##  PTGIS
  ##  COL2A1
  ##  PTX3
  write.table(intersection_module4, file=outputfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
    
    
    
  }else if(length(data3)>1){
  
    ## The data in the inputfile4 and inputfile5 are empty, but the data in the other three inputfiles are not.
    ## So we plot a Venn digaram to compare these three data sets(inputfile1, inputfile2, and inputfile3).
    
    
    ## plot the Venn Diagram
    pdf(outputPlot)
    out.plot <-draw.triple.venn(
      area1 = length(data1),
      area2 = length(data2), 
      area3 = length(data3),
      n12=length(intersect(data1,data2)),
      n23=length(intersect(data2,data3)),
      n13=length(intersect(data1,data3)),
      n123=length(Reduce(intersect, list(data1,data2,data3))),
      category = c("input1","input2","input3"),
      fill = c("blue","red", "green"),
      lty="solid",
      cex =2, 
      cat.cet=2,
      cat.col =c("blue","red", "green")
    )
    print(out.plot)
    dev.off()
    ## export the intersection module detected over the intersection threshold.
    ## get common area of certain input data set in the Venn Digram.
    ## common area between genes from the inputfile1 and inputfile2. 
    node12=as.character(intersect(data1,data2))
    ## common area between genes from the inputfile2 and inputfile3. 
    node23=as.character(intersect(data2,data3))
    ## common area between genes from the inputfile1 and inputfile3. 
    node13=as.character(intersect(data1,data3))
    ## common area among genes from the inputfile1, inputfile2 and inputfile3.
    node123=as.character(Reduce(intersect, list(data1,data2,data3)))
    ## detect the intersection module detected by at least one alogrithm, which is a combination of all nodes detected by these three algorithms.
    module_over_one=unique(c(data1,data2,data3))
    ## detect the intersection module detected by over any two alogrithms, which contains the intersection module detected by any two algortihsm, or all these three algorithms.
    module_over_two=unique(c(node12, node23, node13, node123))
    ## detect the intersection module detected by all three alogrithms
    module_over_three=unique(c(node123))
    ## get the intersection module detected by over certain number of algorithms,  which is defined by the intersectionThreshold.
    intersection_module3=switch(intersectionThreshold,module_over_one,module_over_two,module_over_three)
    
    ##  The output of the intersection module is a list of genes  which are present in the intersection module as follow:
    
    ##  RBP1
    ##  IRF7
    ##  ID4
    ##  PTGIS
    ##  COL2A1
    ##  PTX3
    write.table(intersection_module3, file=outputfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
    
  }else if(length(data2)>1)
 
  { ## The data in the inputfile3, inputfile4 and inputfile5 are empty, but the data in the other two inputfiles are not.
    ## we plot a Venn digaram to compare these two data sets(inputfile1, and inputfile2).
    
    pdf(outputPlot)
    
    out.plot <-  draw.pairwise.venn(
      area1 = length(data1),
      area2 = length(data2),
      cross.area = length(intersect(data1,data2)),
      category = c("input1", "input2"),
      fill = c("orange", "red"),
      lty = "blank",
      cex = 2,
      cat.cex = 2,
      cat.pos = c(0, 0),
      cat.dist = 0.03,
      cat.just = list(c(0, 0), c(0, 0)),
      ext.pos = 30,
      ext.dist = -0.05,
      ext.length = 0.85,
      ext.line.lwd = 2,
      ext.line.lty = "dashed"
    )
    
    print(out.plot)
    dev.off()
    ## detect the intersection module detected by at least one alogrithm, which is a combination of all nodes detected by these two algorithms.
    module_over_one=unique(c(data1,data2))
    ## detect the intersection module detected by all two alogrithms
    module_over_two=unique(c(intersect(data1,data2)))
    ## get the intersection module detected by over certain number of algorithms,  which is defined by the intersectionThreshold.
    intersection_module2=switch(intersectionThreshold,module_over_one,module_over_two)
    ##  The output of the intersection module is a list of genes  which are present in the intersection module as follow:
    
    ##  RBP1
    ##  IRF7
    ##  ID4
    ##  PTGIS
    ##  COL2A1
    ##  PTX3
    
    write.table(intersection_module2, file=outputfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
    
  }else{
    # Only one input was imported and can not perform the comparison
  
    print('Only one data imported and can not perform the comparison!')
  }
}
