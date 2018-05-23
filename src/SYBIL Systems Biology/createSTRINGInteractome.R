##
## Create STRING interactome
##
suppressPackageStartupMessages(library(RGalaxy))

createSTRINGInteractome <- function(
  stringData           = GalaxyInputFile(),
  stringProteinAliases = GalaxyInputFile(),
  confidenceThreashold = GalaxyNumericParam(0.7, required=TRUE, min=0, max=1),
  interactome          = GalaxyOutput("interactome", "tabular")
) {
  ## This function is to retrieve protein-protein interaction from data downloaded from STRING database
  ## There were two input files(protein-protein interaction record and the aliases file)
  ## 1. The protein-protein interaction records, which is stored in the text file 'resources/10090.protein.actions.test.txt'
  ## The format of interaction record were below:
  ## item_id_a	item_id_b	mode	action	a_is_acting	score
  ##  10090.ENSMUSP00000000001	10090.ENSMUSP00000000153	binding		0	540
  ##  10090.ENSMUSP00000000001	10090.ENSMUSP00000000369	activation	activation	1	229
  ##  column 1 is  the STRING id of node A, column 2 is the STRING id of node B, 
  ##  column 3 is the mode of interaction, column 4 is the action nature,
  ##  column 5 ("a_is_acting")is the direction of interaction, if column 5 is 1, which means the interaction is from node A to node B
  ##  if column 5 is 0, it means the interaction is from node B to node A
  ##  column 6 is the confidence score. if the score is more than the confidence trheshold, the interaction are regared as interaction at high confidence score
  
  ## 2. The aliases file maped STRING id to the offical gene symbol
  ## The format of aliase file are:
  ##   ZFP524	10090.ENSMUSP00000083533
  ## column 1 is the offical gene name of the STRING id, 
  ## column 2 is the STRING id
  ## This aliases file currently is manually retrieved from STRING database by following method:
  ## we chose search by multiple names in the web of STRING database,
  ## uploade a text file contained STRING id, in which one STRING id per row,
  ## then search it, and chose to save the format of network in "other formate". STING will promote a new page for the saved result
  ## Nextly, we chose the network description option in the bottom of STRING saved result which were in a text format.
  ## We selected the first two collumns of the network descritpion, which will be stored in the aliases file.
  
  
  ## output needs to be in the format as the createSTRINGInteractome
  ## which  is 2 columns of gene names
  ## Read in the data 
  ppi = read.table(stringData, header=TRUE, sep='\t')

  ## Split into two depending on direction
  ppiA2B = ppi[ppi$a_is_acting==0,]
  ppiB2A = ppi[ppi$a_is_acting==1,]

  ## Rename columns according to the direction
  names(ppiA2B)[names(ppiA2B)=="item_id_a"] <- "from"
  names(ppiA2B)[names(ppiA2B)=="item_id_b"] <- "to"

  names(ppiB2A)[names(ppiB2A)=="item_id_b"] <- "from"
  names(ppiB2A)[names(ppiB2A)=="item_id_a"] <- "to"

  ## merge the results
  ppi = rbind(ppiA2B,ppiB2A)

  ## cut out the columns we require
  ppi <- ppi[c(1,2,3,6)]
  

  ## filter the interactome using the confidence threashold
  ppi<-ppi[ppi$score>=(confidenceThreashold*1000),]
  ppi1<-ppi[c(1,2)]
 
  ## We also need to replace the Gene names using the stringProteinAliases data
  ## the output format is 2 columns of gene names
  aliases=read.table(stringProteinAliases, header=FALSE, sep='\t')
  ## convert all STRING ID in the firt two columns of ppi by gene symbol annotaed in stringProteinAliases
  ##ppi2<-cbind(as.vector( aliases[match(ppi1[,1],aliases[,2]),1]),as.vector( aliases[match(ppi1[,2],aliases[,2]),1]))
  ## create two new columns in pppi1
  ppi1[,3]<-c("NaN")
  ppi1[,4]<-c("NaN")
  ## insert the gene symbol in aliases file whose STRING id matched the STRING id in the column item_id_a of ppi1
  ppi1[,3]<-as.character( aliases[match(ppi1[,1],aliases[,2]),1])
  ## insert the gene symbol in aliases file whose STRING id matched the STRING id in the column item_id_b of ppi1 
  ppi1[,4]<-as.character( aliases[match(ppi1[,2],aliases[,2]),1])
  
  

  ## output should be two columns of gene symbol
  ## remove the first two columns of ppi1
  ppi2<-ppi1[c(3,4)]
  ## remove rows contains NA,
  ppi2<-na.omit(ppi2)
  write.table(ppi2, file=interactome, row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA", sep="\t",append=TRUE)
 
 
 
}