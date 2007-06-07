\name{transform}
\alias{ptransform}

\title{Transform geographic coordinates from one projection into another}
\description{
Transformation of geographics coordinates from one projection to
another, using PROJ.4 library.
}
\usage{
ptransform(data, src.proj, dst.proj, silent=TRUE)
}

\arguments{
  \item{data}{ input (list, matrix or 2d-array) - see details below. }
  \item{src.proj}{ description of the source projection}
  \item{dst.proj}{ description of the destination projection}
  \item{silent}{ if set to \code{TRUE}, warnings will be suppressed }
}
\details{
  The data can be a list of vectors (if the list contains more than
  three entries, only first three entries are used and a warning is
  issued), a two-dimensional matrix or array (the number of columns
  or rows must be at most three). If the input is a list then the result
  will be a list with the entries named \code{x}, \code{y} and \code{z},
  otherwise the result is a matrix with three columns.

  When the list form is used, inputs are recycled with a warning when
  necessary. All unspecified coordinates are set to zero.

  Note that datum files must be installed in order to be able to perform
  datum shifts.
  
  Full details of projection arguments available from website below.
}
\value{
  A matrix with three columns or list with projected coordinates. If the
  input was a list then the output will be a list, otherwise a matrix.
}
\references{\url{http://www.remotesensing.org/proj/}}
\author{Simon Urbanek}

\examples{
data(state)
sc <- cbind(state.center$x, state.center$y)
## this is essentially the same as project except
## that the default lat/long input is in radians
tr <- ptransform(sc/180*pi, '+proj=latlong', '+proj=merc')
## we can compare it with the project result
res <- project(sc, c(proj="merc"))
## ptransform has z coordinate which is 0 for this projection
summary(tr - cbind(res, 0))
}
\keyword{spatial}
