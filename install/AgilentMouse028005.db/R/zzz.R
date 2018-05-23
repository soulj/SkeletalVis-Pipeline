datacache <- new.env(hash=TRUE, parent=emptyenv())

AgilentMouse028005 <- function() showQCData("AgilentMouse028005", datacache)
AgilentMouse028005_dbconn <- function() dbconn(datacache)
AgilentMouse028005_dbfile <- function() dbfile(datacache)
AgilentMouse028005_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
AgilentMouse028005_dbInfo <- function() dbInfo(datacache)

AgilentMouse028005ORGANISM <- "Mus musculus"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "AgilentMouse028005.sqlite", package=pkgname, lib.loc=libname)
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
    ann_objs <- createAnnObjs.SchemaChoice("MOUSECHIP_DB", "AgilentMouse028005", "chip AgilentMouse028005", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("AgilentMouse028005.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(AgilentMouse028005_dbconn())
}

