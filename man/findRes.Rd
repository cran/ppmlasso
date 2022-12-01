\name{findRes}
\alias{findRes}
\title{Choose spatial resolution for analysis}
\description{This function produces a plot to choose the optimal spatial resolution for
analysis. A point process model is calculated for each nominated spatial resolution and the
log-likelihood of all fitted models are plotted against the spatial resolutions.}
\usage{
findRes(scales, lambda = 0, coord = c("X", "Y"), sp.xy, env.grid, 
formula, tol = 0.01, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scales}{A vector of spatial resolutions for which to produce the plot.}
  \item{lambda}{The penalty for each fitted spatial resolution. This should be a
single value such that only one point process model is computed for each spatial resolution.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates,
used by \code{\link{getEnvVar}}.}
  \item{sp.xy}{A matrix of species locations containing at least one column representing
longitude and one column representing latitude, as in \code{\link{ppmlasso}}.}
  \item{env.grid}{The geo-referenced matrix of environmental grids, as in \code{\link{ppmlasso}}.}
  \item{formula}{The formula of the fitted model, as in \code{\link{ppmlasso}}.}
  \item{tol}{An optional argument to specify the tolerance level of coordinate error passed to an internal call to the \code{\link{griddify}} function, set to 0.01 by default.}
  \item{\dots}{Further arguments passed to \code{\link{ppmlasso}}.
}
}
\details{This function produces a plot which can be used to judge an optimal spatial resolution for analysis. As the spatial resolution gets finer, the log-likelihood tends to stabilise to a constant value. The largest spatial resolution at which the log-likelihood appears to stabilise may be considered optimal for model fitting.}
\value{A plot of log-likelihood versus spatial resolution.}
\references{Renner, I.W. et al (2015). Point process models for presence-only analysis. \emph{Methods in Ecology and Evolution} \bold{6}, 366-379.}
\author{Ian W. Renner}
\examples{
data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
scales = c(0.5, 1, 2, 4, 8, 16)
ppm.form = ~ poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2)
findRes(scales, formula = ppm.form, sp.xy = sub.euc, env.grid = sub.env)
}