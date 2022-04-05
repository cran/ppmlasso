\name{sampleQuad}
\alias{sampleQuad}
\title{Generate regular grid of quadrature points with environmental data}
\description{This function generates a regular grid of quadrature points and
associated environmental data at a nominated spatial resolution.}
\usage{
sampleQuad(env.grid, sp.scale, coord = c("X", "Y"), quadfilename = "Quad", tol = 0.01, 
writefile = TRUE)
}
\arguments{
  \item{env.grid}{The geo-referenced matrix of environmental grids. It must
have a vector of longitude and a vector of latitude.}
  \item{sp.scale}{The spatial resolution at which to sample quadrature points.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates.}
  \item{quadfilename}{An optional argument containing the prefix of the name of the saved file. The
default is \code{"Quad"} so that a matrix generated at a spatial resolution of 1 would
be saved in the file \code{"Quad1.RData"}. A file is saved for every resolution given
in \code{sp.scale}.}
  \item{tol}{An optional argument to specify the tolerance level of coordinate error passed to an internal call to the \code{\link{griddify}} function, set to 0.01 by default.}
  \item{writefile}{A logical argument to determine whether the output should be written to a file or not. If \code{TRUE} (the default), the output will be saved with the file name as described in the \code{quadfilename} argument.}
}
\value{The output is a matrix of quadrature points at the spatial resolution
supplied to \code{sp.scale}. If a vector of resolutions is supplied, the output is a list of
matrices for each spatial resolution.}
\author{Ian W. Renner}
\examples{
data(BlueMountains)
quad.1 = sampleQuad(env.grid = BlueMountains$env, sp.scale = 1, quadfilename = NA, 
writefile = FALSE)
}
\keyword{Quadrature}