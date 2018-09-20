suppressPackageStartupMessages(library(RGalaxy))

TF_Enrichment <- function(
  differentialExpression = GalaxyInputFile(required=TRUE),
  species = GalaxyCharacterParam(c("Human","Mouse")),
  foldChangeOnly = GalaxyLogicalParam(),
  foldchange = GalaxyNumericParam(1.5),
  padj = GalaxyNumericParam(0.05),
  enrichmentTable = GalaxyOutput("TF_EnrichmentTable", "tabular")) {
  
  
  suppressPackageStartupMessages(library("RcisTarget"))
  suppressPackageStartupMessages(library("dplyr"))
  
  #load the data
  differentialExpression<-read.delim(differentialExpression)
  
  colnames(differentialExpression)[1]<-"GeneSymbol"
  
  if (foldChangeOnly==TRUE){
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      slice(which.max(abs(2))) %>% as.data.frame
    
    
    differentiallyExpressedGenes<-na.omit(differentialExpression[ abs(differentialExpression[,2])>=log2(foldchange),])
    geneList<-list(geneSet=differentiallyExpressedGenes[,1])
    
  } else{
    
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      slice(which.min(3)) %>% as.data.frame
    differentiallyExpressedGenes<-na.omit(differentialExpression[ abs(differentialExpression[,2])>=log2(foldchange) & differentialExpression[,3] <=padj,])
    geneList<-list(geneSet=differentiallyExpressedGenes[,1])
  }
  
  #run the TF enrichment analysis
  if (species=="Human"){
    
    suppressPackageStartupMessages(library("RcisTarget.hg19.motifDatabases.20k"))
    data(hg19_10kbpAroundTss_motifRanking)
    data(hg19_direct_motifAnnotation)
    data(hg19_inferred_motifAnnotation)
    motifRankings <- hg19_10kbpAroundTss_motifRanking
    geneList$geneSet <- geneList$geneSet [geneList$geneSet %in% getRanking(motifRankings)$rn]
    motifs_AUC <- calcAUC(geneList, motifRankings, nCores=1)
    motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3, 
                                               motifAnnot_direct=hg19_direct_motifAnnotation,
                                               motifAnnot_inferred=hg19_inferred_motifAnnotation)
    motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                       rankings=motifRankings, 
                                                       geneSets=geneList)
    
  } else {
    
    suppressPackageStartupMessages(library("RcisTarget.mm9.motifDatabases.20k"))
    data(mm9_10kbpAroundTss_motifRanking)
    data(mm9_direct_motifAnnotation)
    data(mm9_inferred_motifAnnotation)
    motifRankings <- mm9_10kbpAroundTss_motifRanking
   # geneList$geneSet <- geneList$geneSet [geneList$geneSet %in% getRanking(motifRankings)$rn]
    motifs_AUC <- calcAUC(geneList, motifRankings, nCores=1)
    motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3, 
                                               motifAnnot_direct=mm9_direct_motifAnnotation,
                                               motifAnnot_inferred=mm9_inferred_motifAnnotation)
    motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                       rankings=motifRankings, 
                                                       geneSets=geneList)
    }
  
  write.table(motifEnrichmentTable_wGenes,file=enrichmentTable,col.names=T,row.names=F,sep="\t",quote=F)

}

