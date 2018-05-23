
\name{ArrayXHuman_dbconn}

\alias{ArrayXHuman_dbconn}
\alias{ArrayXHuman_dbfile}
\alias{ArrayXHuman_dbschema}
\alias{ArrayXHuman_dbInfo}


\title{Collect information about the package annotation DB}

\description{
  Some convenience functions for getting a connection object to (or collecting
  information about) the package annotation DB.
}

\usage{
  ArrayXHuman_dbconn()
  ArrayXHuman_dbfile()
  ArrayXHuman_dbschema(file="", show.indices=FALSE)
  ArrayXHuman_dbInfo()
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
  \code{ArrayXHuman_dbconn} returns a connection object to the
  package annotation DB.  IMPORTANT: Don't call
  \code{\link[DBI:dbConnect]{dbDisconnect}} on the connection object
  returned by \code{ArrayXHuman_dbconn} or you will break all the
  \code{\link[AnnotationDbi:AnnDbObj-class]{AnnDbObj}} objects defined
  in this package!

  \code{ArrayXHuman_dbfile} returns the path (character string) to the
  package annotation DB (this is an SQLite file).

  \code{ArrayXHuman_dbschema} prints the schema definition of the
  package annotation DB.

  \code{ArrayXHuman_dbInfo} prints other information about the package
  annotation DB.
}

\value{
  \code{ArrayXHuman_dbconn}: a DBIConnection object representing an
  open connection to the package annotation DB.

  \code{ArrayXHuman_dbfile}: a character string with the path to the
  package annotation DB.

  \code{ArrayXHuman_dbschema}: none (invisible \code{NULL}).

  \code{ArrayXHuman_dbInfo}: none (invisible \code{NULL}).
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
  dbGetQuery(ArrayXHuman_dbconn(), "SELECT COUNT(*) FROM probes")

  ArrayXHuman_dbschema()

  ArrayXHuman_dbInfo()
}

\keyword{utilities}
\keyword{datasets}

