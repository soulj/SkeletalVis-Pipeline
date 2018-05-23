\name{AgilentMouse014868ORGANISM}
\alias{AgilentMouse014868ORGANISM}
\alias{AgilentMouse014868ORGPKG}

\title{The Organism information for AgilentMouse014868}
\description{
  AgilentMouse014868ORGANISM is an R object that contains a single item: a
  character string that names the organism for which AgilentMouse014868 was
  built.
  AgilentMouse014868ORGPKG is an R object that contains a chararcter vector
  with the name of the organism package that a chip package depends on
  for its gene-centric annotation.
}
\details{
  Although the package name is suggestive of the organism for which it
  was built, AgilentMouse014868ORGANISM provides a simple way to programmatically
  extract the organism name.
  AgilentMouse014868ORGPKG provides a simple way to programmatically
  extract the name of the parent organism package.  The parent organism
  package is a strict dependency for chip packages as this is where the
  gene cetric information is ultimately extracted from.  The full
  package name will always be this string plus the extension ".db".  But
  most programatic acces will not require this extension, so its more
  convenient to leave it out.
}

\seealso{
  \itemize{
    \item \code{\link[AnnotationDbi]{AnnotationDb-class}} for use of
          the \code{select()} interface.
  }
}

\examples{
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.

## Bimap interface:
AgilentMouse014868ORGANISM
AgilentMouse014868ORGPKG
}

\keyword{datasets}

