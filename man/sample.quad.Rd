\name{sample.quad}
\alias{sample.quad}
\title{Generate regular grid of quadrature points with environmental data}
\description{This function generates a regular grid of quadrature points and
associated environmental data at a nominated spatial resolution.}
\usage{
sample.quad(env.grid, sp.scale, coord = c("X", "Y"), file = "Quad")
}
\arguments{
  \item{env.grid}{The geo-referenced matrix of environmental grids. It must
have a vector of longitude and a vector of latitude.}
  \item{sp.scale}{The spatial resolution at which to sample quadrature points.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates.}
  \item{file}{An optional argument containing the prefix of the name of the saved file. The
default is \code{"Quad"} so that a matrix generated at a spatial resolution of 1 would
be saved in the file \code{"Quad1.RData"}. A file is saved for every resolution given
in \code{sp.scale}.}
}
\value{The output is a matrix of quadrature points at the spatial resolution
supplied to \code{sp.scale}. If a vector of resolutions is supplied, the output is a list of
file names containing the saved matrices of quadrature points stored as \code{dat.quad}.}
\author{Ian W. Renner}
\examples{
data(BlueMountains)
quad.1 = sample.quad(env.grid = BlueMountains$env, sp.scale = 1, file = "Quad")
}
\keyword{Quadrature}
