\name{AgilentRat028279CHR}
\alias{AgilentRat028279CHR}
\title{Map Manufacturer IDs to Chromosomes}
\description{
 AgilentRat028279CHR is an R object that provides mappings between a
 manufacturer identifier and the chromosome that contains the gene of interest. 
}

\details{
  Each manufacturer identifier maps to a vector of chromosomes. Due to
  inconsistencies that may exist at the time the object was built,
  the vector may contain more than one chromosome (e.g., the identifier may map
  to more than one chromosome). If the chromosomal location is unknown,
  the vector will contain an \code{NA}.
  
  Mappings were based on data provided by: Entrez Gene  
  ftp://ftp.ncbi.nlm.nih.gov/gene/DATA  
  With a date stamp from the source of: 2016-Sep26
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
x <- AgilentRat028279CHR
# Get the probe identifiers that are mapped to a chromosome
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  # Get the CHR for the first five probes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
}
\keyword{datasets}

