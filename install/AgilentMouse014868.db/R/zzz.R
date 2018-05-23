datacache <- new.env(hash=TRUE, parent=emptyenv())

AgilentMouse014868 <- function() showQCData("AgilentMouse014868", datacache)
AgilentMouse014868_dbconn <- function() dbconn(datacache)
AgilentMouse014868_dbfile <- function() dbfile(datacache)
AgilentMouse014868_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
AgilentMouse014868_dbInfo <- function() dbInfo(datacache)

AgilentMouse014868ORGANISM <- "Mus musculus"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "AgilentMouse014868.sqlite", package=pkgname, lib.loc=libname)
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
    ann_objs <- createAnnObjs.SchemaChoice("MOUSECHIP_DB", "AgilentMouse014868", "chip AgilentMouse014868", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("AgilentMouse014868.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(AgilentMouse014868_dbconn())
}

