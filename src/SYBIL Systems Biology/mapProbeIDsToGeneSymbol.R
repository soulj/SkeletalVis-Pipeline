#
# Lookup the gene symbol for a probe ID, using the DAVID web service.
#
# Uses biomart:https://www.bioconductor.org/packages/3.3/bioc/html/biomaRt.html
# Source: https://github.com/Bioconductor-mirror/DAVIDQuery/

## This function is to map expressed probe id to their official gene symbol.



suppressPackageStartupMessages(library(RGalaxy))

mapProbeIDsToGeneSymbol <- function(
  inputfile1       = GalaxyInputFile(required=TRUE),
  
  probeColumnIndex = GalaxyNumericParam(1, required=TRUE, min=1),
  outputfile1      = GalaxyOutput("geneMappedData", "tabular")
) {
  ## read the expressed probe id from the inputfile1, 'resources/probe_id_test_data.txt'.
  data <- read.table(inputfile1, sep='\t', as.is=TRUE)
 
 
  ## The input of data type is various, for example, it can be alist of probe ids (the example 1,data from the file, 'resources/probe_id_test_data.txt') as follow:
  
  ##  1415670_at
  ##  1415671_at
  ##  1415672_at
  ##  1415673_at
  ##  1415674_a_at
  
  ## or it could be a table of which one column shows the probe ids (the example 2, data from the file, 'resources/MLIIAnalysisSample.txt')
 
  ##   1457666_s_at	-7.02139374728657	6.58070004035722	-63.3535154086283	9.36195572844208e-10	2.26583618442668e-05	12.5909129914816
  ##   1453084_s_at	-6.44840680955418	7.85291297507557	-59.9592657471613	1.30478859328421e-09	2.26583618442668e-05	12.3924174957527
  ##   1448690_at	-5.18688284519683	8.71616517265489	-51.7315400072555	3.17588359738665e-09	2.26583618442668e-05	11.8085964105249
  ##   1440962_at	-4.78008804033848	9.18120867658225	-46.6543450157693	5.91721193757442e-09	2.26583618442668e-05	11.3564704075447
  ##   1456073_s_at	-7.94914679211289	8.65262733810269	-46.6267762205736	5.93831379672011e-09	2.26583618442668e-05	11.3537840808768
  ##   1457589_at	-4.83882176893325	8.89268286346451	-45.0403467422846	7.31468078513742e-09	2.26583618442668e-05	11.19454684087
  
  
  ## probeColumnIndex is 1 in the example 1 and the example 2.The first one column contains all expressed probes.
  
  presentProbe<-data[,probeColumnIndex]
  
  probeid<-as.character(presentProbe)
  ## remove duplicate probe ids from the list.
  presentProbe1<-unique(probeid)
  
  ## use biomaRt to map present probe ids to official gene symbol.
  library(biomaRt)
  ## connect the selected database, "ENSEMBL_MART_ENSEMBL"and data set,"mmusculus_gene_ensembl".
  ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host="www.ensembl.org")
  ## retireve infromation from the biomaRt database
  ## The offical gene symbol is in the column, "mgi_symbol".
  biomaRtData<-getBM(attributes=c('affy_mouse430_2','mgi_symbol'),filters = 'affy_mouse430_2', values = presentProbe1, mart = ensembl)
  
  
  ## The retrieved the biomaRtData  of the example 1 is a table as follow:
  ##    affy_mouse430_2   mgi_symbol
  ## 1       1460301_at      Vmn1r62
  ## 2       1434928_at       Gas2l1
  ## 3       1450270_at      Pcdhb11
  ## 4       1456710_at      Pcdhb11
  ## 5     1437423_a_at         Sra1
  ## 6       1424457_at         Sra1
  
  ## The retrieved the biomaRtData  of theexample 2 is a table as follow:
  
  ##    affy_mouse430_2    mgi_symbol
  ##  1      1460301_at       Vmn1r62
  ##  2      1422182_at         Hnf4g
  ##  3      1450518_at         Hnf4g
  ##  4      1434928_at        Gas2l1
  ##  5      1441422_at        Gas2l1
  ##  6      1430595_at        Gas2l1
  
 
  
  ## remove redundant rows that different probe ids mapped to the same gene symbol.
  biomaRtData2<-subset(biomaRtData, !duplicated(biomaRtData[,"mgi_symbol"]))
  
  
  ## The table of biomaRtData2 for the example 1 should be as follow:
    
  ##    affy_mouse430_2     mgi_symbol
  ## 1       1460301_at        Vmn1r62
  ## 2       1434928_at         Gas2l1
  ## 3       1450270_at        Pcdhb11
  ## 4     1437423_a_at           Sra1
  ## 5     1428154_s_at       Ppapdc1b
  ## 6       1417572_at            Mpg
  
  ## The table of biomaRtData2 for the example 2 should be as follow:
  ##     affy_mouse430_2    mgi_symbol
  ##  1       1460301_at       Vmn1r62
  ##  2       1422182_at         Hnf4g
  ##  4       1434928_at        Gas2l1
  ##  7       1422617_at        Gm2092
  ##  9       1429854_at 2900092C05Rik
  ##  10      1418680_at      Serpind1
  
  ## replace the probe ids with the mapped gene symbol according to biomaRtData2 and convert all gene symbols to its uppercase.
 
  data[,probeColumnIndex]<-toupper(as.character( biomaRtData2[match(data[,probeColumnIndex],biomaRtData2[,"affy_mouse430_2"]),"mgi_symbol"]))
  ## remove rows whose gene symbol column is "NA".
  data2<-data[!is.na(data[,probeColumnIndex]),]
  
  ## the output of the example 1 should be as follow:
  ##         V1
  ## 1    COPG1
  ## 2 ATP6V0D1
  ## 3   GOLGA7
  ## 4     PSPH
  ## 5  TRAPPC4
  ## 6     DPM2
  
  
  ## the output of the example 2 should be as follow:
  ##        V1        V2       V3        V4           V5           V6       V7
  ## 1 COL22A1 -6.448407 7.852913 -59.95927 1.304789e-09 2.265836e-05 12.39242
  ## 2   KCNK1 -5.186883 8.716165 -51.73154 3.175884e-09 2.265836e-05 11.80860
  ## 3  SLC8A3 -4.780088 9.181209 -46.65435 5.917212e-09 2.265836e-05 11.35647
  ## 4   PANX3 -7.949147 8.652627 -46.62678 5.938314e-09 2.265836e-05 11.35378
  ## 5    FAT3 -4.838822 8.892683 -45.04035 7.314681e-09 2.265836e-05 11.19455
  ## 6 COL24A1 -4.547539 6.566030 -43.43194 9.104990e-09 2.265836e-05 11.02328
  
  write.table(data2, file=outputfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, na="", sep="\t")
  
}
