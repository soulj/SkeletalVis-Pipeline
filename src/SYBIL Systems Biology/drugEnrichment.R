suppressPackageStartupMessages(library("RGalaxy"))

drugEnrichment <- function(
  differentialExpression = GalaxyInputFile(required=TRUE),
  homology = GalaxyInputFile(required=FALSE),
  foldChangeOnly         = GalaxyLogicalParam(),
  species=GalaxySelectParam(c("Human","Mouse","Rat","Horse","Zebrafish","Cow","Pig")),
  padj = GalaxyNumericParam(0.05),
  enrichedDrugsMimic    = GalaxyOutput("enrichedDrugsMimic", "tabular"),
  enrichedDrugsReverse    = GalaxyOutput("enrichedDrugsReverse", "tabular")) {
  
  suppressPackageStartupMessages(library("httr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library("org.Mm.eg.db"))
  suppressPackageStartupMessages(library("org.Rn.eg.db"))
  suppressPackageStartupMessages(library("org.Dr.eg.db"))
  suppressPackageStartupMessages(library("org.Ss.eg.db"))
  suppressPackageStartupMessages(library("org.Bt.eg.db"))
  suppressPackageStartupMessages(library("org.Ecaballus.eg.db"))

  
  #Get the enriched drugs and their targets from a list of differentially expressed genes
  getEnrichedDrugs<-function(diffExp,species,homology){
    
    #map to human orthologs if necessary
    if (species != "Human") {
      homology<-read.delim(homology,stringsAsFactors = F)
      diffExp<-diffExp[ diffExp[,1] %in% homology[,3],]
      homology<-homology[ homology[,3] %in% diffExp[,1],-1:-2]
      colnames(homology)[2]<-"gene_name.human"
      diffExp <- merge(diffExp,homology,by.x=1,by.y=1)
      diffExp[,1]<-diffExp[,"gene_name.human"]
    }

    upGenes <- diffExp[ diffExp[,2]>=log2(1.5),1]
    downGenes <- diffExp[ diffExp[,2]<=log2(1/1.5),1]
    
    
    url <- "http://amp.pharm.mssm.edu/L1000CDS2/query"
    
    genes <- list("upGenes"=upGenes,"dnGenes"=downGenes)
    config.mimic <- list("aggravate"=TRUE,"searchMethod"="geneSet","share"=FALSE,"combination"=FALSE,"db-version"="latest")
    config.reverse <- list("aggravate"=FALSE,"searchMethod"="geneSet","share"=FALSE,"combination"=FALSE,"db-version"="latest")
    payload.mimic <- list("data"=genes,"config"=config.mimic)
    payload.reverse <- list("data"=genes,"config"=config.reverse)
    response.mimic <- POST(url, body = payload.mimic, encode = "json")
    response.reverse <- POST(url, body = payload.reverse, encode = "json")
    drugs.mimic<-parseResponse(response.mimic)
    drugs.reverse<-parseResponse(response.reverse)
    
    #get the targets of the enriched drugs
    drugs.mimic$targets<-sapply(drugs.mimic$pubchem_id,getDrugTargets)
    drugs.reverse$targets<-sapply(drugs.reverse$pubchem_id,getDrugTargets)
    
    drugs.mimic$pert_desc<-apply(drugs.mimic,1,function(x) {
      if(as.character(x[1])=="-666"& !is.na(x[7])){
        getDrugName(as.character(x[7]))
      } else{
        return(as.character(x[1]))
      }
    })
    
    drugs.reverse$pert_desc<-apply(drugs.reverse,1,function(x) {
      if(as.character(x[1])=="-666" & !is.na(x[7])){
        getDrugName(as.character(x[7]))
      } else{
        return(as.character(x[1]))
      }
    })
      

    
    return(list(drugs.mimic=drugs.mimic,drugs.reverse=drugs.reverse))
  }
  
  #parse the response from the L1000CDS2 API
  parseResponse<-function(response){
    
    response <- httr::content(response)[[2]]
    
    #extract the overlapping up and down genes - convert to a string
    overlap.up<-sapply(response,function(x) paste(as.character(unlist(x[["overlap"]][[1]])), sep="' '", collapse=", "))
    overlap.down<-sapply(response,function(x) paste(as.character(unlist(x[["overlap"]][[2]])), sep="' '", collapse=", "))
    
    #create a table of the results
    response<-lapply(response,"[",-10)
    response<-bind_rows(lapply(response,function(y){as.data.frame(t(unlist(y)),stringsAsFactors=FALSE)}))
    response$pert_dose<-paste0(response$pert_dose,response$pert_dose_unit)
    response$pert_time<-paste0(response$pert_time,response$pert_time_unit)
    response<-response[,c("pert_desc","pert_id","cell_id","pert_dose","pert_time","score","pubchem_id")]
    response$overlap.up<-overlap.up
    response$overlap.down<-overlap.down
    
    return(response)
  }
  
  getDrugName<-function(CID){
    url<-paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",CID,"/description/JSON")
    r<-GET(url)
    r<-httr::content(r)
    
    return(r[[1]][[1]][[1]][["Title"]])
  }

  
  #get the targets of the drugs by CID
  getDrugTargets<-function(CID){
    
    url<-paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",CID,"/assaysummary/CSV")
    r<-GET(url)
    
    if(status_code(r)!=200) {
      return(NA)
    }
    
    r<-httr::content(r)
    #get the gene symbols of the active outcomes
    r<-na.omit(r[ r$`Bioactivity Outcome`=="Active","Target GeneID"])
    
    #convert the entrezID to the human/mouse gene symbol
    if(length(r$`Target GeneID`)>0){
      targets<-c()
      if (any(r$`Target GeneID` %in% keys(org.Hs.eg.db))){
        targets<-na.omit(AnnotationDbi::select(org.Hs.eg.db,as.character(unique(r$`Target GeneID`)),c("SYMBOL")))
      }
      if (any(r$`Target GeneID` %in% keys(org.Mm.eg.db))){
        targets<-rbind(targets,na.omit(AnnotationDbi::select(org.Mm.eg.db,as.character(unique(r$`Target GeneID`)),c("SYMBOL"))))
      }
      targets<-sort(unique(targets$SYMBOL))
      targets<-paste(targets, sep="' '", collapse=", ")
      return(targets)
    }else{
      return(NA)
    }
    
    
  }
  
  
  differentialExpression<-read.delim(differentialExpression,stringsAsFactors = F)
  
  colnames(differentialExpression)[1]<-"GeneSymbol"
  
  if (foldChangeOnly==TRUE){
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.max(abs(2))) %>% as.data.frame
    differentialExpression.sig<-na.omit(differentialExpression[ abs(differentialExpression[,2])>=log2(1.5),])
  } else{
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.min(3)) %>% as.data.frame
    differentialExpression.sig<-na.omit(differentialExpression[ abs(differentialExpression[,2])>=log2(1.5) & differentialExpression[,3] <=padj,])
  }
  
  #get the enriched pathways
  enrichedDrugsTable<- getEnrichedDrugs(differentialExpression.sig,species,homology)
 
  write.table(enrichedDrugsTable[[1]],file=enrichedDrugsMimic,col.names=T,row.names=F,sep="\t",quote=F)
  write.table(enrichedDrugsTable[[2]],file=enrichedDrugsReverse,col.names=T,row.names=F,sep="\t",quote=F)
  
}
