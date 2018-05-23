datacache <- new.env(hash=TRUE, parent=emptyenv())

AgilentHuman039494 <- function() showQCData("AgilentHuman039494", datacache)
AgilentHuman039494_dbconn <- function() dbconn(datacache)
AgilentHuman039494_dbfile <- function() dbfile(datacache)
AgilentHuman039494_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
AgilentHuman039494_dbInfo <- function() dbInfo(datacache)

AgilentHuman039494ORGANISM <- "Homo sapiens"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "AgilentHuman039494.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    txdb <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"ChipDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, txdb, envir=ns)
    namespaceExport(ns, dbNewname)
        
    ## Create the AnnObj instances
    ann_objs <- createAnnObjs.SchemaChoice("HUMANCHIP_DB", "AgilentHuman039494", "chip AgilentHuman039494", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("AgilentHuman039494.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(AgilentHuman039494_dbconn())
}

