#!/usr/bin/env Rscript

## begin warning handler
withCallingHandlers({

library(methods) # Because Rscript does not always do this

options('useFancyQuotes' = FALSE)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RGalaxy"))


option_list <- list()

option_list$accessionNumber <- make_option('--accessionNumber', type='character')
option_list$sampleTable <- make_option('--sampleTable', type='character')
option_list$downloadedFiles <- make_option('--downloadedFiles', type='character')


opt <- parse_args(OptionParser(option_list=option_list))




getRNASeqExpressionData <- function (accessionNumber = GalaxyCharacterParam(), sampleTable = GalaxyInputFile(required = T, 
    formatFilter = "tabular"), downloadedFiles = GalaxyOutput("files", 
    "tabular")) 
{
    library(ENAbrowseR)
    library(curl)
    sampleTable <- read.delim(sampleTable)
    results <- ena_search(query = paste0("study_accession=", 
        accessionNumber), result = "read_run")
    results <- results[results$fastq_ftp != "", ]
    fastqs <- unlist(strsplit(results$fastq_ftp, ";"))
    fastqs <- unique(grep(paste(sampleTable$File, collapse = "|"), 
        fastqs, value = TRUE))
    fastqs <- paste0("ftp://", fastqs)
    fileNames <- sapply(fastqs, basename)
    dir.create("fastqFiles")
    sapply(fastqs, function(fastq) curl_fetch_disk(fastq, paste0("fastqFiles/", 
        basename(fastq))))
    write.table(fileNames, file = downloadedFiles, col.names = T, 
        row.names = F, sep = "\t", quote = F)
}

params <- list()
for(param in names(opt))
{
    if (!param == "help")
        params[param] <- opt[param]
}

setClass("GalaxyRemoteError", contains="character")
wrappedFunction <- function(f)
{
    tryCatch(do.call(f, params),
        error=function(e) new("GalaxyRemoteError", conditionMessage(e)))
}


suppressPackageStartupMessages(library(RGalaxy))
do.call(getRNASeqExpressionData, params)

## end warning handler
}, warning = function(w) {
    cat(paste("Warning:", conditionMessage(w), "\n"))
    invokeRestart("muffleWarning")
})
