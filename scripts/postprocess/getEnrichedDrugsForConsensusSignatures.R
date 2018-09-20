getEnrichedDrugs<-function(upGenes,downGenes){
  
  print("hello")
  
  suppressPackageStartupMessages(library("httr"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library("org.Mm.eg.db"))
  
  url <- "http://amp.pharm.mssm.edu/L1000CDS2/query"
  
  if(is.null(upGenes)|is.null(downGenes)){
    return(list(NA,NA))
  }
  
  genes <- list("upGenes"=upGenes,"dnGenes"=downGenes)
  config.mimic <- list("aggravate"=TRUE,"searchMethod"="geneSet","share"=FALSE,"combination"=FALSE,"db-version"="latest")
  config.reverse <- list("aggravate"=FALSE,"searchMethod"="geneSet","share"=FALSE,"combination"=FALSE,"db-version"="latest")
  payload.mimic <- list("data"=genes,"config"=config.mimic)
  payload.reverse <- list("data"=genes,"config"=config.reverse)
  response.mimic <- POST(url, body = payload.mimic, encode = "json")
  response.reverse <- POST(url, body = payload.reverse, encode = "json")
  drugs.mimic<-parseResponse(response.mimic)
  drugs.reverse<-parseResponse(response.reverse)
  
  if (is.null(drugs.mimic)|is.null(drugs.reverse)) return(list(NA,NA))
  
  #get the targets of the enriched drugs
  # drugs.mimic$targets<-sapply(drugs.mimic$pubchem_id,getDrugTargets)
  # drugs.reverse$targets<-sapply(drugs.reverse$pubchem_id,getDrugTargets)
  
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
  
  response <- try(httr::content(response)[[2]])
  if(class(response) == "try-error") return(NULL)
  
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