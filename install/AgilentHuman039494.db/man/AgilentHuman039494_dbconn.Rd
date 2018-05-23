
\name{AgilentHuman039494_dbconn}

\alias{AgilentHuman039494_dbconn}
\alias{AgilentHuman039494_dbfile}
\alias{AgilentHuman039494_dbschema}
\alias{AgilentHuman039494_dbInfo}


\title{Collect information about the package annotation DB}

\description{
  Some convenience functions for getting a connection object to (or collecting
  information about) the package annotation DB.
}

\usage{
  AgilentHuman039494_dbconn()
  AgilentHuman039494_dbfile()
  AgilentHuman039494_dbschema(file="", show.indices=FALSE)
  AgilentHuman039494_dbInfo()
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
  \code{AgilentHuman039494_dbconn} returns a connection object to the
  package annotation DB.  IMPORTANT: Don't call
  \code{\link[DBI:dbConnect]{dbDisconnect}} on the connection object
  returned by \code{AgilentHuman039494_dbconn} or you will break all the
  \code{\link[AnnotationDbi:AnnDbObj-class]{AnnDbObj}} objects defined
  in this package!

  \code{AgilentHuman039494_dbfile} returns the path (character string) to the
  package annotation DB (this is an SQLite file).

  \code{AgilentHuman039494_dbschema} prints the schema definition of the
  package annotation DB.

  \code{AgilentHuman039494_dbInfo} prints other information about the package
  annotation DB.
}

\value{
  \code{AgilentHuman039494_dbconn}: a DBIConnection object representing an
  open connection to the package annotation DB.

  \code{AgilentHuman039494_dbfile}: a character string with the path to the
  package annotation DB.

  \code{AgilentHuman039494_dbschema}: none (invisible \code{NULL}).

  \code{AgilentHuman039494_dbInfo}: none (invisible \code{NULL}).
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
  dbGetQuery(AgilentHuman039494_dbconn(), "SELECT COUNT(*) FROM probes")

  AgilentHuman039494_dbschema()

  AgilentHuman039494_dbInfo()
}

\keyword{utilities}
\keyword{datasets}

