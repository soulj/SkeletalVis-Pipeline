datacache <- new.env(hash=TRUE, parent=emptyenv())

ArrayXHuman <- function() showQCData("ArrayXHuman", datacache)
ArrayXHuman_dbconn <- function() dbconn(datacache)
ArrayXHuman_dbfile <- function() dbfile(datacache)
ArrayXHuman_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
ArrayXHuman_dbInfo <- function() dbInfo(datacache)

ArrayXHumanORGANISM <- "Homo sapiens"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "ArrayXHuman.sqlite", package=pkgname, lib.loc=libname)
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
    ann_objs <- createAnnObjs.SchemaChoice("HUMANCHIP_DB", "ArrayXHuman", "chip ArrayXHuman", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("ArrayXHuman.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(ArrayXHuman_dbconn())
}

