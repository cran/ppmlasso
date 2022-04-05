\name{getEnvVar}
\alias{getEnvVar}
\title{Extract environmental data to presence locations}
\description{Given a matrix of quadrature points and a list of species presences, this
function extracts environmental data to presence locations using bilinear interpolation.}
\usage{
getEnvVar(sp.xy, env.grid, env.scale, coord = c("X", "Y"), envfilename = "SpEnvData", 
tol = 0.01, writefile = TRUE)
}
\arguments{
  \item{sp.xy}{A matrix of species locations containing at least one column representing
longitude and one column representing latitude.}
  \item{env.grid}{The geo-referenced matrix of environmental grids.}
  \item{env.scale}{The spatial resolution of the quadrature scheme from which the 
environmental data is extracted.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates,
as in \code{\link{sampleQuad}}.}
  \item{envfilename}{An optional argument containing the name of the saved file. Setting
\code{envfilename =} \code{"SpEnvData"} will save a matrix \code{sp.dat} containing the species presence
locations and the interpolated environmental data to the file \code{"SpEnvData.RData"}.}
  \item{tol}{An optional argument to specify the tolerance level of coordinate error passed to an internal call to the \code{\link{griddify}} function, set to 0.01 by default.}
  \item{writefile}{A logical argument to determine whether the output should be written to a file or not. If \code{TRUE} (the default), the output will be saved with the file name as described in the \code{envfilename} argument.}
}
\details{At a given species location with coordinates \eqn{(x, y)}, the interpolated value of the
environmental variable \eqn{z} is calculated as a weighted average of \eqn{z} at four reference
quadrature points \eqn{(x^{(1)}, y^{(1)})}, \eqn{(x^{(1)}, y^{(2)})}, \eqn{(x^{(2)}, y^{(1)})} and \eqn{(x^{(2)}, y^{(2)})}
that form a square of nominated side length \code{env.scale} surrounding \eqn{(x, y)}.}
\value{A matrix containing locations of species presences in the first two columns
and the interpolated environmental data in the remaining columns.}
\author{Ian W. Renner}
\examples{
data(BlueMountains)
species.env = getEnvVar(BlueMountains$eucalypt, env.grid = BlueMountains$env, env.scale = 0.5,
envfilename = NA, writefile = FALSE)
}
\keyword{Interpolation}