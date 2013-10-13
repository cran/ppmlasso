\name{findres}
\alias{findres}
\title{Choose spatial resolution for analysis}
\description{This function produces a plot to choose the optimal spatial resolution for
analysis. A point process model is calculated for each nominated spatial resolution and the
log-likelihood of all fitted models are plotted against the spatial resolutions.}
\usage{
findres(scales, lambda = 0, coord = c("X", "Y"), sp.xy, env.grid, formula, ...)
}
\arguments{
  \item{scales}{A vector of spatial resolutions for which to produce the plot.}
  \item{lambda}{The LASSO penalty for each fitted spatial resolution. This should be a
single value such that only one point process model is computed for each spatial resolution.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates,
used by \code{\link{env.var}}.}
  \item{sp.xy}{A matrix of species locations containing at least one column representing
longitude and one column representing latitude, as in \code{\link{ppmlasso}}.}
  \item{env.grid}{The geo-referenced matrix of environmental grids, as in \code{\link{ppmlasso}}.}
  \item{formula}{The formula of the fitted model, as in \code{\link{ppmlasso}}.}
  \item{\dots}{Further arguments passed to \code{\link{ppmlasso}}.}
}
\value{A plot of log-likelihood versus spatial resolution.}
\author{Ian W. Renner}
\examples{
data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
scales = c(0.5, 1, 2, 4, 8, 16)
ppm.form = ~ poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree=2)
findres(scales, formula = ppm.form, sp.xy = sub.euc, env.grid = sub.env)
}
