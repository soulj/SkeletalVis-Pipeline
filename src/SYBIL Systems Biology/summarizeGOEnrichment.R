#
# Uses the bioconductor library to perform GO enrichment analysis for significant genes.
#
#
# See https://www.bioconductor.org/packages/release/bioc/manuals/affy/man/affy.pdf
# CEL format: http://www.affymetrix.com/estore/support/developer/powertools/changelog/gcos-agcc/cel.html.affx
#

suppressPackageStartupMessages(library(RGalaxy))

summarizeGOEnrichment <- function(
         topicGeneInputfile = GalaxyInputFile(required=TRUE),
  backgroundSampleInputfile = GalaxyInputFile(required=TRUE),
            goTermTypeIndex = GalaxyNumericParam(1, required=TRUE, min=1),
          goTermOutputfile1 = GalaxyOutput("output", "tabular")
) {
  suppressPackageStartupMessages(library(affy))
  suppressPackageStartupMessages(library(biomaRt))
  suppressPackageStartupMessages(library(GO.db))
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  suppressPackageStartupMessages(library(topGO))
  suppressPackageStartupMessages(library(clusterProfiler))

  ## This function will perform GO enrichment analysis for the topic genes we are interested in.
  ## There are three inputs required for the function:topic gene input file,background sample input file and gene ontology term type index.
  ## The topic gene input file is a text file which lists the topic genes we are interested in.
  ## It is the target gene set for the GO enrichment analysis.
  ## The format of topicGeneInputfile contained one gene symbol per row as below:
  
  ## COL2A1
  ## CXCL10
  ## F3
  ## FMOD
  ## IBSP
  ## ID4
  
  ## The background sample input file is a CEL files. It includes all genes annotated to a GO term in the entire background set. 
  ## The GO term type index is a selection of the GO term type. 
  ## It has three values: "1" indicates the biological process, "2" indicates the cellular component and "3" indicates the molecular function.
  
  ## read the topic gene input file.
  data <- read.table(topicGeneInputfile, sep='\t', as.is=TRUE)
  gene_symbol<-as.character(data[,1])
  
  # ReadAffy supports a character vector of files to read into the batch.
  ## The backgroundSampleInputfile includeded  one  samples below:
  ## 'resources/GO/GSM1072318_pes1_KI_OB_I.CEL', which is alread stored in the SYBIL database.
  print('Reading affy batch..')
  affyInput <- c(backgroundSampleInputfile)
  affyInput <- affyInput[affyInput != 'None']
  print(affyInput)
  dataAffyBatch <- ReadAffy(filenames=affyInput)
  ## extract all probe ids of the sample as the background gene set.
  id1<-as.character(probeNames(dataAffyBatch))
  ## The gene symbols were converted to the entrezgene id for GO enrichment analysis by the biomaRt.
  library(biomaRt)
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
  ## The topic genes were mapped to their entrezgene id.
  data2<-getBM(attributes=c('mgi_symbol','entrezgene'),filters = 'mgi_symbol', values = gene_symbol, mart = ensembl)
  topic_gene<-na.omit(data2$entrezgene)
  ## The probe ids of the background gene set were mapped to their entrezgene id.
  data3<-getBM(attributes=c('affy_mouse430_2', 'entrezgene','mgi_symbol'),filters = 'affy_mouse430_2', values = id1, mart = ensembl)
  all_gene<-na.omit(data3$entrezgene)
  ## if goTermTypeIndex is 1, the GO terms of biological process will be analyzed.
  ## if goTermTypeIndex is 2, the GO terms of cellular component will be analyzed.
  ## if goTermTypeIndex is 2, the GO terms of molecular function will be analyzed.
  if(goTermTypeIndex=="1")
         {goTerms <- enrichGO(gene= as.character(topic_gene),"org.Mm.eg.db",ont = "BP",pAdjustMethod = "BH", universe=as.character(all_gene),pvalueCutoff  = 0.1,minGSSize  = 5,readable =FALSE)}
  else if(goTermTypeIndex=="2")      
  {goTerms <- enrichGO(gene= as.character(topic_gene),"org.Mm.eg.db",ont = "CC",pAdjustMethod = "BH", universe=as.character(all_gene),pvalueCutoff  = 0.1,minGSSize  = 5,readable =FALSE)}
  else if(goTermTypeIndex=="3")      
   {goTerms <- enrichGO(gene= as.character(topic_gene),"org.Mm.eg.db",ont = "MF",pAdjustMethod = "BH", universe=as.character(all_gene),pvalueCutoff  = 0.1,minGSSize  = 5,readable =FALSE)}
  ## Calculate the statistic results of the GO term found.
  output<-summary(goTerms)
  
  write.table(output, file=goTermOutputfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, na="", sep="\t", append=TRUE)
  
  ## The output will contain nie columns: ID of GO term, Description of GO term,	GeneRatio, BgRatio,	pvalue, p.adjust,	qvalue, 	geneID and	Count.
  ## It will be a statistic table as follow:
  
  ## ID	|Description |	GeneRatio	|BgRatio|	pvalue |	p.adjust |	qvalue |	geneID |	Count
  ## GO:0031988	membrane-bounded vesicle	85/221	3028/17556	4.07056946905975e-14	4.12579790709624e-12	3.12211726173667e-12	103711/105348/108156/108664/108679/11364/11502/11671/116847/11769/11972/11993/12183/12263/12462/12468/12469/12527/12847/12876/13194/14066/14229/14376/15891/16592/16783/17523/17758/17975/18570/18571/18607/18791/18807/18826/19156/19173/19184/19349/20103/20133/20224/20501/20615/20778/20983/21787/22022/22027/22152/22166/22240/22329/22334/22629/227753/23918/24070/246703/26433/26440/27054/27370/433375/52502/54130/54161/56421/56434/56443/56451/56463/57437/60409/64294/65114/66092/66111/66201/66687/66848/68135/68539/71665	85
  ## GO:0070062	extracellular vesicular exosome	76/221	2515/17556	4.09900940716639e-14	4.12579790709624e-12	3.12211726173667e-12	103711/105348/108156/108664/108679/11364/11502/11671/116847/11769/11972/11993/12183/12263/12462/12468/12469/12527/12847/12876/13194/14066/14229/14376/16592/16783/17523/17758/17975/18570/18571/18791/18807/18826/19156/19173/19184/19349/20103/20133/20224/20501/20778/21787/22022/22027/22152/22166/22329/22334/22629/227753/23918/24070/246703/26433/26440/27370/433375/52502/54130/56421/56434/56443/56451/56463/57437/64294/65114/66092/66201/66687/66848/68135/68539/71665	76
  ## GO:1903561	extracellular vesicle	76/221	2515/17556	4.09900940716639e-14	4.12579790709624e-12	3.12211726173667e-12	103711/105348/108156/108664/108679/11364/11502/11671/116847/11769/11972/11993/12183/12263/12462/12468/12469/12527/12847/12876/13194/14066/14229/14376/16592/16783/17523/17758/17975/18570/18571/18791/18807/18826/19156/19173/19184/19349/20103/20133/20224/20501/20778/21787/22022/22027/22152/22166/22329/22334/22629/227753/23918/24070/246703/26433/26440/27370/433375/52502/54130/56421/56434/56443/56451/56463/57437/64294/65114/66092/66201/66687/66848/68135/68539/71665	76
  
  ## geneID is the entrezgene id
  
}
