suppressPackageStartupMessages(library("RGalaxy"))

GigaSubnetworks <- function(
  differentialExpression = GalaxyInputFile(required=TRUE),
  interactome            = GalaxyInputFile(required=TRUE),
  foldChangeOnly         = GalaxyLogicalParam(),
  species=GalaxySelectParam(c("Human","Mouse","Rat","Horse","Zebrafish","Cow","Pig")),
  moduleNodes            = GalaxyOutput("moduleNodes", "text"),
  modulePlots            = GalaxyOutput("modulePlots", "pdf"),
  visNetworks = GalaxyOutput("visNetworks", "rdata"),
  summaryTable = GalaxyOutput("summaryTable", "text")
  
) {
  
  suppressPackageStartupMessages(library("igraph"))
  suppressPackageStartupMessages(library("squash"))
  #devtools::install_github("briatte/ggnet")
  suppressPackageStartupMessages(library("ggnet"))
  suppressPackageStartupMessages(library("intergraph"))
  suppressPackageStartupMessages(library("visNetwork"))
  suppressPackageStartupMessages(library("goseq"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library("org.Mm.eg.db"))
  suppressPackageStartupMessages(library("org.Rn.eg.db"))
  suppressPackageStartupMessages(library("org.Dr.eg.db"))
  suppressPackageStartupMessages(library("org.Ss.eg.db"))
  suppressPackageStartupMessages(library("org.Bt.eg.db"))
  suppressPackageStartupMessages(library("org.Ecaballus.eg.db"))
  
  getTopGOTerms<-function(geneLists,geneLengths,differentialExpression,species){
    
    #map from gene symbols to REACTOME
    if (species == "Human") {
      gene2GO<-AnnotationDbi::select(org.Hs.eg.db,keys(org.Hs.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Mouse"){
      gene2GO<-AnnotationDbi::select(org.Mm.eg.db,keys(org.Mm.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Cow"){
      gene2GO<-AnnotationDbi::select(org.Bt.eg.db,keys(org.Bt.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Horse"){
      gene2GO<-AnnotationDbi::select(org.Ecaballus.eg.db,keys(org.Ecaballus.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Zebrafish"){
      gene2GO<-AnnotationDbi::select(org.Dr.eg.db,keys(org.Dr.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else if (species == "Pig"){
      gene2GO<-AnnotationDbi::select(org.Ss.eg.db,keys(org.Ss.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    } else{
      gene2GO<-AnnotationDbi::select(org.Rn.eg.db,keys(org.Rn.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL")
    }
    gene2GO<-unstack(gene2GO[,c(1,5)])
    
    TopGOTerms<-sapply(geneLists,getTopGOTerm,geneLengths,differentialExpression,gene2GO)
    return(TopGOTerms)
  }
  
  getTopGOTerm<-function(geneList,geneLengths,differentialExpression,gene2GO){
    
    genes <- ifelse(differentialExpression[,1] %in% geneList ,1,0)
    genes <- as.data.frame(genes)
    colnames(genes) <- "DEgenes"
    genes$bias.data <- 1
    genes$pwf <- 1
    rownames(genes) <- differentialExpression[,1]

    goTerms <- goseq(genes,gene2cat=gene2GO,method="Hypergeometric")
    
    return(goTerms[1,"term"])
    
  }
  
  
  getVisSubnetworks<-function(subnetTable,nodeList,network){
    
    #set up the colours
    maxVal<-max(abs(subnetTable[,3])) + 0.05
    minVal<--maxVal
    
    breaksList = seq(minVal, maxVal, by = 0.05)
    colfunc <- colorRampPalette(c("green", "white","red"))
    map <- makecmap(subnetTable[,3],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
    
    subnetworks<-lapply(nodeList,function(x) induced.subgraph(network,x))
    
    visNetworks<-lapply(subnetworks,getVisSubnetwork,map,subnetTable)
    
    return(visNetworks)
    
  }
  
  getVisSubnetwork<-function(subnetwork,colours,subnetTable){
    
    FCs<-subnetTable[match(V(subnetwork)$name,subnetTable[,1]),3]
    colours<-cmap(FCs, map = colours)
    
    nodes<-data.frame(id=V(subnetwork)$name,title=as.character(FCs),label=V(subnetwork)$name,color=colours,size=16,font.size=22)
    edges<-as.data.frame(get.edgelist(subnetwork))
    edges<-data.frame(from=edges[,1],to=edges[,2],color = "black")
    
    network<-visNetwork(nodes,edges)
    return(network)
  }
  
  
  plotSubnetworks<-function(subnetTable,nodeList,network){
    
    #set up the colours
    maxVal<-max(abs(subnetTable[,3])) + 0.05
    minVal<--maxVal
    
    print(maxVal)
    print(minVal)
    breaksList = seq(minVal, maxVal, by = 0.05)
    colfunc <- colorRampPalette(c("green", "white","red"))
    map <- makecmap(subnetTable[,3],n = 3,breaks = breaksList,symm = T,colFn = colfunc)
    
    subnetworks<-lapply(nodeList,function(x) induced.subgraph(network,x))
    subnetworks<-lapply(subnetworks,simplify)
    plots<-lapply(subnetworks,plotSubnetwork,map,subnetTable)
    
    return(plots)
    
  }
  
  plotSubnetwork<-function(subnetwork,colours,subnetTable){
    FCs<-subnetTable[match(V(subnetwork)$name,subnetTable[,1]),3]
    colours<-cmap(FCs, map = colours)
    layout <- layout.fruchterman.reingold(subnetwork)
    if(length(layout[,1])==2){
      #layout <- layout.random(subnetwork)
    }
    g<-ggnet2(subnetwork,mode=layout,color=colours,label = TRUE,size=6,edge.color = "black",edge.size = 0.3,label.size = 4, layout.exp = 0.2)+
      theme(plot.margin=unit(c(1,1,1,1),"cm"))
    return(g)
  }


  #R implementation of the the GIGA perl script - couple orders of magnitude faster, but ugly.
  #finds sub-networks from a ranked list
  runGIGA=function(expressionData,inputGraph,max_number=20) {
    Scores<-data.frame(name=expressionData[,1],geneRank=rank(-expressionData$Pi),1:nrow(expressionData))
    adj<-get.adjacency(inputGraph,sparse=T)
    neighboursAdj<-lapply(1:nrow(adj),getNeighboursAdj,adj)
    neighboursAdjSelf<-mapply(c,1:nrow(adj),neighboursAdj)
    
    Scores<-Scores[order(Scores[,3]),]
    #get the seeds
    seedCluster<-initaliseCluster(Scores,inputGraph,neighboursAdj,neighboursAdjSelf)
    #expand the from the seeds
    expandedCluster<-clusterExpand(Scores,seedCluster,completedCluster=c(),inputGraph,max_number,neighboursAdjSelf)
    
    #keep adding nodes until max size or pvalue doesn't improve
    while(length(expandedCluster[["cluster"]])>0) {
      seedCluster<-addNewMin(Scores,expandedCluster,inputGraph,neighboursAdjSelf)
      expandedCluster<-clusterExpand(Scores,seedCluster,seedCluster[["completedCluster"]],inputGraph,max_number,neighboursAdjSelf)
      
    }
    #return the final sub-networks
    finalClusters<-getFinalClusters(Scores,expandedCluster,inputGraph)
    clusters<-finalClusters[[1]]
    pvalues<-finalClusters[[2]]
    clusters<-lapply(clusters,function(x) V(inputGraph)$name[x])
    return(list(clusters=clusters,pvalues=pvalues))
    
  }
  
  
  getNeighboursAdj<-function(i,adj){
    neighbours<-which(adj[i,]>0)
    return(neighbours)
    
  }
  
  
  #function to find the local min in the PPI network
  initaliseCluster = function(Scores,inputGraph,neighboursAdj,neighboursAdjSelf) {
    
    
    #identify the neihbours of each node and find the local minima - those with no connections to a node with a higher rank.
    size<-length(neighboursAdjSelf)
    NodeDataRankNamed<-Scores[,3]
    names(NodeDataRankNamed)=Scores$geneRank
    neighbours<-lapply(neighboursAdjSelf,function(x) NodeDataRankNamed[x])
    localMinIndex<-lapply(neighbours, function(x) which.min(as.numeric(names(x))))
    localMin<-which(localMinIndex==1)
    
    #expand the initial cluster
    
    cluster<-as.list(localMin)
    names(cluster)<-V(inputGraph)$name[localMin]
    oldcluster<-cluster
    cluster<-lapply(cluster,function(x) NodeDataRankNamed[x])
    oldpvalue<-rep(1,length(cluster))
    names(oldpvalue)<-names(cluster)
    clustermax<-localMin
    
    neighbours<-neighboursAdj[localMin]
    neighbours<-sapply(neighbours,function(x) NodeDataRankNamed[x])
    neighbours<-sapply(neighbours,function(x) x[sort.list(as.numeric(names(x)))])
    neighbours<-lapply(neighbours,'[',1)
    cluster<-mapply(c, cluster, neighbours, SIMPLIFY=FALSE)
    clustermax<-sapply(cluster,function(x) names(x[2]))
    
    return(list(clustermax=clustermax,oldpvalue=oldpvalue,cluster=cluster,oldcluster=oldcluster))
    
  }
  
  
  #function to expand the existing sub-network
  clusterExpand = function(Scores,seedCluster,completedCluster,inputGraph,max_number,neighboursAdjSelf) {
    NodeDataRankNamed=Scores[,3]
    names(NodeDataRankNamed)<-Scores$geneRank
    clustermax<-seedCluster[["clustermax"]]
    oldpvalue<-seedCluster[["oldpvalue"]]
    cluster<-seedCluster[["cluster"]]
    oldcluster<-seedCluster[["oldcluster"]]
    
    cluster_new<-cluster
    tobePValued<-c()
    
    size<-vcount(inputGraph)
    
    #keep expanding while you can (neighbours less than current max rank)
    while (length(cluster_new)>0) {
      
      neigbours<-getNeigbours(cluster_new,neighboursAdjSelf)
      lastexpansion<-cluster_new
      
      #assign the expanded cluster with the rank
      neigbours<-lapply(neigbours,function(x) NodeDataRankNamed[x])
      #extract from 1 to the max rank in each cluster
      clustermax<-as.numeric(clustermax)
      expansion<-lapply(seq_along(neigbours),function(x)  neigbours[[x]][(as.numeric(names( neigbours[[x]]))<=clustermax[x])])
      names(expansion)<-names(cluster_new)
      #remove any clusters that are too big
      expansion_tooBig<-lapply(expansion,function(x) length(x)<max_number)
      cluster_new<-expansion[unlist(expansion_tooBig)]
      completedCluster<-c(completedCluster,oldcluster[!unlist(expansion_tooBig)])
      oldcluster<-oldcluster[unlist( expansion_tooBig)]
      clustermax<-clustermax[unlist(expansion_tooBig)]
      #seperate the clusters that have changed in size
      clusterlength<-lapply(cluster_new,length)
      oldclusterlength<-lapply(lastexpansion,length)
      continueCluster<-lapply(names(cluster_new),function(x) clusterlength[[x]]>oldclusterlength[[x]])
      if (length(continueCluster>1)) {
        tobePValued<-c(tobePValued,cluster_new[!unlist(continueCluster)])
      }
      cluster_new<-cluster_new[unlist(continueCluster)]
      clustermax<-clustermax[unlist(continueCluster)]
      oldcluster<-oldcluster[unlist(continueCluster)]
    }
    
    if (length(tobePValued)>0) {
      oldcluster<-seedCluster[["oldcluster"]]  
      oldcluster<-oldcluster[match(names(tobePValued),names(oldcluster))]
      oldpvalue<-oldpvalue[names(oldpvalue) %in% names(tobePValued)]
      currentPValue<-lapply(tobePValued,getPvalue,size)
      expansion_improved<-lapply(names(currentPValue),function(x) currentPValue[[x]]<oldpvalue[[x]])
      names(expansion_improved)<-names(currentPValue)
      cluster<-tobePValued[unlist(expansion_improved)]
      worseClusters<-oldcluster[!unlist(expansion_improved)]
      completedCluster<-c(completedCluster,worseClusters)
      oldpvalue<-currentPValue[unlist(expansion_improved)]

    }
    else {cluster<-NULL}
    
    return(list(oldpvalue=oldpvalue,cluster=cluster,completedCluster=completedCluster))
  }
  
  
  addNewMin= function(Scores,expandedCluster,inputGraph,neighboursAdjSelf) {
    
    NodeDataRankNamed<-Scores[,3]
    names(NodeDataRankNamed)<-Scores$geneRank
    cluster<-expandedCluster[["cluster"]]
    oldpvalue<-expandedCluster[["oldpvalue"]]
    completedCluster<-expandedCluster[["completedCluster"]]
    
    oldcluster<-cluster
    newNeighbours<-getNeigbours(cluster,neighboursAdjSelf)
    
    
    clustermax<-c()
    newMin<-lapply(newNeighbours,function(x) NodeDataRankNamed[x])
    names(newMin)<-names(cluster)
    newMin<-lapply(newMin,function(x) x[sort.list(as.numeric(names(x)))])
    addition<-lapply(seq_along(newMin),function(x) newMin[[x]][!(newMin[[x]] %in% cluster[[x]])][1]) 
    cluster<-mapply(c, cluster, addition, SIMPLIFY=FALSE)
    clustermax<-sapply(cluster,function(x) max(as.numeric(names(x))))
    
    return(list(clustermax=clustermax,oldpvalue=oldpvalue,cluster=cluster,oldcluster=oldcluster,completedCluster=completedCluster))
    
  }
  
  getPvalue=function(x,size) {
    x<-as.numeric(names(x))
    maxRank<-max(x)
    pvalue<-maxRank/size
    if (length(x)==1) return(pvalue)
    for (k in 1:(length(x)-1)) {
      pvalue<-pvalue*((maxRank-k)/(size-k))
    }
    return(pvalue)
  }
  
  getNeigbours = function(cluster,neighboursAdj) {
    
    newNeighbours<-lapply(cluster,function(x) unlist(neighboursAdj[ unlist(x)]))
    newNeighbours<-lapply(newNeighbours,function(x) x[!duplicated(x)])
    
    return(newNeighbours)
  }
  
  getFinalClusters=function(Scores,expandedCluster,inputGraph){
    size<-vcount(inputGraph)  
    NodeDataRankNamed<-Scores[,3]
    names(NodeDataRankNamed)<-Scores$geneRank
    completedCluster<-expandedCluster[["completedCluster"]]
    completedCluster<-lapply(completedCluster,function(x) NodeDataRankNamed[x])
    completedClusterpValue<-lapply(completedCluster, getPvalue,size)
    completedCluster<-completedCluster[completedClusterpValue<1/size]
    completedClusterpValue<-lapply(completedCluster, getPvalue,size)
    completedClusterpValue<-data.frame(name=names(completedCluster),pvalue=unlist(completedClusterpValue))
    completedClusterpValue<-completedClusterpValue[order(completedClusterpValue$pvalue),]
    completedClusterpValue$localmin<-match(completedClusterpValue$name,V(inputGraph)$name)
    
    localmin<-completedClusterpValue[,3]
    finalCluster<-c()
    localminused<-c()
    alreadydone<-c()
    finalPvalues<-c()
    
    
    
    for (i in 1:nrow(completedClusterpValue)) {
      if ((completedClusterpValue[i,3] %in% localminused)==FALSE) {
        candidateCluster<-completedCluster[[as.character(completedClusterpValue[i,1])]]
        localminpresent<-candidateCluster[candidateCluster %in% localmin ]
        localminused<-c(localminpresent,localminused)
        finalCluster<-c(finalCluster,list(candidateCluster))
        finalPvalues<-c(finalPvalues,completedClusterpValue[i,"pvalue"])
      }
    }
    
    
    return(list(finalCluster,finalPvalues))
    
  }
  
  #load the data
  differentialExpression<-read.delim(differentialExpression)
  
  interactome <- read.table(interactome, sep='\t', as.is=TRUE)
  interactome <- graph.data.frame(interactome, directed=F)
  
 # geneLengths<-read.delim(geneLengths)
  
  #calculate the Pi value for scoring the nodes as approriate
  colnames(differentialExpression)[1]<-"GeneSymbol"
  differentialExpression$absFC<-abs(differentialExpression[,2])

  if(foldChangeOnly == T) {
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.max(abs(2))) %>% as.data.frame
    differentialExpression$Pi<-differentialExpression$absFC
    
  } else {
    differentialExpression<-  differentialExpression %>%
      group_by(GeneSymbol) %>%
      dplyr::slice(which.min(3)) %>% as.data.frame
    differentialExpression$logAdjPval<- -log10(differentialExpression[,3])
    differentialExpression$Pi<-(differentialExpression$absFC*differentialExpression$logAdjPval)
    
  }
  
  #filter the network based on expressed genes and extract the largest connected component
  differentialExpression<-na.omit(differentialExpression)
  presentList<-na.omit(match(differentialExpression[,1],V(interactome)$name))
  interactome<-induced.subgraph(interactome,presentList)
  interactome<-decompose.graph(interactome)
  interactome<-interactome[[which.max(sapply(interactome, vcount))]]
  
  #filter and order the expression table based on genes in the network
  presentList<-na.omit(match(V(interactome)$name,differentialExpression[,1]))
  differentialExpression<-differentialExpression[presentList,]
  
  #run GIGA
  results<-runGIGA(expressionData = differentialExpression,inputGraph = interactome,max_number = 20)
  pvalues<-results$pvalues

  
  #output the results
  output <- data.frame(subnetworks = rep(1:length(results$clusters), sapply(results$clusters, length)), gene_name = unlist(results$clusters))
  output<-merge(output, differentialExpression, by.x = 'gene_name', by.y = 1, all.x=T)
  output<-output[order(output$subnetworks, output$gene_name),]
  write.table(output[, !(names(output) %in% c('absFC', 'logAdjPval', 'Pi'))], moduleNodes, sep = '\t', row.names = F,quote = F)
  
  #summaryTable
  topGOTerms<-getTopGOTerms(results$clusters,geneLengths,differentialExpression,species)
  resultsSummary<-data.frame(Network=1:length(results$clusters),Size=sapply(results$clusters,length),pvalue=results$pvalues,topGOBPTerm=topGOTerms)
  write.table(resultsSummary, summaryTable, sep = '\t', row.names = F,quote = F)
  
  #interactive networks for shiny app
  networks<-getVisSubnetworks(output,results$clusters,interactome)
  save(networks,file=visNetworks)
  
  # Plot the subnetworks
  suppressPackageStartupMessages(library("ggpubr"))
  plots<-plotSubnetworks(output,results$clusters,interactome)
  pdf(file=modulePlots)
  arrange<-ggarrange(plotlist=plots,ncol=1,nrow=2)
  lapply(arrange,plot)
  dev.off()
  

}
