
\name{AgilentRat014879_dbconn}

\alias{AgilentRat014879_dbconn}
\alias{AgilentRat014879_dbfile}
\alias{AgilentRat014879_dbschema}
\alias{AgilentRat014879_dbInfo}


\title{Collect information about the package annotation DB}

\description{
  Some convenience functions for getting a connection object to (or collecting
  information about) the package annotation DB.
}

\usage{
  AgilentRat014879_dbconn()
  AgilentRat014879_dbfile()
  AgilentRat014879_dbschema(file="", show.indices=FALSE)
  AgilentRat014879_dbInfo()
}

\arguments{
  \item{file}{
    A connection, or a character string naming the file to print to (see
    the \code{file} argument of the \code{\link[base]{cat}} function for
    the details).
  }
  \item{show.indices}{
    The CREATE INDEX statements are not shown by default.
    Use \code{show.indices=TRUE} to get them.
  }
}

\details{
  \code{AgilentRat014879_dbconn} returns a connection object to the
  package annotation DB.  IMPORTANT: Don't call
  \code{\link[DBI:dbConnect]{dbDisconnect}} on the connection object
  returned by \code{AgilentRat014879_dbconn} or you will break all the
  \code{\link[AnnotationDbi:AnnDbObj-class]{AnnDbObj}} objects defined
  in this package!

  \code{AgilentRat014879_dbfile} returns the path (character string) to the
  package annotation DB (this is an SQLite file).

  \code{AgilentRat014879_dbschema} prints the schema definition of the
  package annotation DB.

  \code{AgilentRat014879_dbInfo} prints other information about the package
  annotation DB.
}

\value{
  \code{AgilentRat014879_dbconn}: a DBIConnection object representing an
  open connection to the package annotation DB.

  \code{AgilentRat014879_dbfile}: a character string with the path to the
  package annotation DB.

  \code{AgilentRat014879_dbschema}: none (invisible \code{NULL}).

  \code{AgilentRat014879_dbInfo}: none (invisible \code{NULL}).
}

\seealso{
  \code{\link[DBI:dbSendQuery]{dbGetQuery}},
  \code{\link[DBI:dbConnect]{dbConnect}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbconn}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbfile}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbschema}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbInfo}}
}

\examples{
  library(DBI)
  ## Count the number of rows in the "probes" table:
  dbGetQuery(AgilentRat014879_dbconn(), "SELECT COUNT(*) FROM probes")

  AgilentRat014879_dbschema()

  AgilentRat014879_dbInfo()
}

\keyword{utilities}
\keyword{datasets}

