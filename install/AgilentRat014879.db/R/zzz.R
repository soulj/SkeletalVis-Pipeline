datacache <- new.env(hash=TRUE, parent=emptyenv())

AgilentRat014879 <- function() showQCData("AgilentRat014879", datacache)
AgilentRat014879_dbconn <- function() dbconn(datacache)
AgilentRat014879_dbfile <- function() dbfile(datacache)
AgilentRat014879_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
AgilentRat014879_dbInfo <- function() dbInfo(datacache)

AgilentRat014879ORGANISM <- "Rattus norvegicus"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "AgilentRat014879.sqlite", package=pkgname, lib.loc=libname)
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
    ann_objs <- createAnnObjs.SchemaChoice("RATCHIP_DB", "AgilentRat014879", "chip AgilentRat014879", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("AgilentRat014879.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(AgilentRat014879_dbconn())
}

