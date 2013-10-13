\name{point.interactions}
\alias{point.interactions}
\title{Calculate point interactions for area-interaction models}
\description{This function calculates point interactions at presence locations and
quadrature points required for fitting a regularisation path of area-interaction
models.}
\usage{
point.interactions(dat.ppm, r, availability = NA)
}
\arguments{
  \item{dat.ppm}{A design matrix generated using the \code{\link{ppm.dat}} function.}
  \item{r}{The radius of point interactions.}
  \item{availability}{An optional binary matrix used in calculating point interactions indicating
whether locations are available (\code{availability = 1}) or not (\code{availability = 0}). 
If no such matrix is provided, \code{availability} is automatically generated with all
values set to \code{1} at a special resolution of half of \code{r}.}
}
\details{Theoretically, the point interaction \eqn{t(y)} at a point \eqn{y} is calculated
as the proportion of available area in a circular region \eqn{Y} or radius \eqn{r} centred
at \eqn{y} that overlaps with circles of radius \eqn{r} centred at other presence locations
(Baddeley & Turner, 2005).

This function discretises the study region at the same spatial resolution
as \code{availability} by defining the matrix \code{occupied}, a fine grid of locations
spanning the study region initialised to zero. The values of \code{occupied} within a 
distance of \code{r} of each presence location are then augmented by \code{1}, such that
\code{occupied} then contains the total number of presence locations with which each
grid location interacts. To prevent unavailable areas from being included in the
calculation of point interactions, the values of \code{occupied} at grid locations for which
\code{availability = 0} are set to zero.

\eqn{t(y)} is then estimated as the proportion of available grid locations within \eqn{Y}
that overlap circular regions around other presence locations.

The availability matrix is particularly useful for regions that have inaccessible areas
(\strong{e.g.} due to the presence of ocean or urban areas).

Finer resolutions of the \code{availability} matrix will yield more precise estimates but
at a cost of greater computation time.}
\value{A vector of point interactions corresponding to the locations contained in the
\code{dat.ppm} argument.}
\references{
Baddeley, A.J. & Turner, R. (2005). Spatstat: an \code{R} package for analyzing spatial
point patterns. \emph{Journal of Statistical Software} \bold{12}, 1-42.
}
\author{Ian W. Renner}
\seealso{
\code{\link{ppmlasso}} for fitting a regularisation path of area-interaction models
}
\examples{
data(BlueMountains)
species.ppm = ppm.dat(sp.xy = BlueMountains$eucalypt, back.xy = BlueMountains$env, 
sp.scale = 1, file.name = "Sp PPM Data") # generate design matrix
species.int = point.interactions(species.ppm, 2, BlueMountains$availability)
}
\keyword{Area-interaction models}
