\name{env.var}
\alias{env.var}
\title{Extract environmental data to presence locations}
\description{Given a matrix of quadrature points and a list of species presences, this
function extracts environmental data to presence locations using bilinear interpolation.}
\usage{
env.var(sp.xy, env.grid, env.scale, coord = c("X", "Y"), file.name = NA)
}
\arguments{
  \item{sp.xy}{A matrix of species locations containing at least one column representing
longitude and one column representing latitude.}
  \item{env.grid}{The geo-referenced matrix of environmental grids.}
  \item{env.scale}{The spatial resolution of the quadrature scheme from which the 
environmental data is extracted.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates,
as in \code{\link{sample.quad}}.}
  \item{file.name}{An optional argument containing the name of the saved file. Setting
\code{file.name =} \code{"Sp Env"} will save a matrix \code{sp.dat} containing the species presence
locations and the interpolated environmental data to the file \code{"Sp Env.RData"}.}
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
species.env = env.var(BlueMountains$eucalypt, env.grid = BlueMountains$env, env.scale = 0.5,
file.name = "Sp Env Data")
}
\keyword{Interpolation}
