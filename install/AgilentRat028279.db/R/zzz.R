datacache <- new.env(hash=TRUE, parent=emptyenv())

AgilentRat028279 <- function() showQCData("AgilentRat028279", datacache)
AgilentRat028279_dbconn <- function() dbconn(datacache)
AgilentRat028279_dbfile <- function() dbfile(datacache)
AgilentRat028279_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
AgilentRat028279_dbInfo <- function() dbInfo(datacache)

AgilentRat028279ORGANISM <- "Rattus norvegicus"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "AgilentRat028279.sqlite", package=pkgname, lib.loc=libname)
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
    ann_objs <- createAnnObjs.SchemaChoice("RATCHIP_DB", "AgilentRat028279", "chip AgilentRat028279", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("AgilentRat028279.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(AgilentRat028279_dbconn())
}

