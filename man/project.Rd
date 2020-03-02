\name{project}
\alias{project}

\title{Projection of coordinates}
\description{
Projection of lat/long coordinates or its inverse.
}
\usage{
project(xy, proj, inverse = FALSE, degrees = TRUE, silent = FALSE,
        ellps.default="sphere")
}

\arguments{
  \item{xy}{ input (list, matrix or 2d-array) - see details below. }
  \item{proj}{ projection definition }
  \item{inverse}{ if \code{TRUE} inverse projection is performed (from a
    cartographic projection into lat/long), otherwise projects from
    lat/long into a cartographic projection. }
  \item{degrees}{ if \code{TRUE} then the lat/long data is assumed to be
    in degrees, otherwise in radians }
  \item{silent}{ if set to \code{TRUE}, warnings will be suppressed }
  \item{ellps.default}{default ellipsoid that will be added if no datum
  or ellipsoid parameter is specified in \code{proj}. Older versions of
  PROJ.4 didn't require a datum (and used sphere by default), but 4.5.0
  and higher always require a datum or an ellipsoid. Set to \code{NA} if
  no datum should be added to \code{proj} (e.g. if you specify an
  ellipsoid directly).}
}
\details{
  The input can be a list of two or more vectors (if the list contains
  more than two entries, only first two entries are used and a warning
  is issued), a two-dimensional matrix or array (the number of columns
  or rows must be exactly two) or a vector of the length 2. For a 2x2
  input the columns are taken as \code{x} and \code{y}. If the input
  is a list then the result will be a list with the entries named
  \code{x} and \code{y}, otherwise the result is a matrix with two
  columns.

  When the list form is used, inputs are recycled with a warning when
  necessary.
  
  \code{proj} specifies the target (or source if inverse)
  projection. The format can be either a single (unnamed) string which
  contains all parameters:
  \code{"+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96"}
  or an unnamed vector of complete individual parameters:
  \code{c("+proj=lcc","+lat_1=33","+lat_2=45","+lat_0=39","+lon_0=-96")}
  or a named vector or list that will be composed into parameters:
  \code{list(proj="lcc", lat_1=33, lat_2=45, lat_0=39, lon_0=-96)}

  if \code{degrees} is \code{TRUE} then the latitude and longitude are
  expected to be in degrees, if \code{FALSE} then in radians.
}
\value{
  A two column matrix or list of coordinates. If the input was a list
  then the output will be a list, otherwise a matrix.
}
\references{\url{https://proj.org/}}

\examples{
## this is just very simple, because we don't want to depend on
## maps package, so we can't show more useful stuff..
data(state)
s <- project(state.center, "+proj=merc")
plot(s, type='n', asp=1)
text(s,, state.abb)
}
\keyword{spatial}
