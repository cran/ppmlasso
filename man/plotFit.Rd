\name{plotFit}
\alias{plotFit}
\title{Plot the predicted intensity of a fitted \code{ppmlasso} model}
\description{This function produces a levelplot of the predicted intensity from a fitted \code{ppmlasso} model.}
\usage{
plotFit(fit, pred.data = data.frame(X = fit$x[fit$pres == 0], 
Y = fit$y[fit$pres == 0], 1, scale(fit$data[fit$pres == 0, -1], 
center = -fit$s.means/fit$s.sds, scale = 1/fit$s.sds)), 
coord = c("X", "Y"), asp = "iso", ylab = "", xlab = "", 
col.regions = heat.colors(1024)[900:1], cuts = length(col.regions), 
cex = 1.4, main.text = paste(toupper(fit$criterion), "fit"), cex.color = 1.4)
}
\arguments{
  \item{fit}{A fitted \code{ppmlasso} object.}
  \item{pred.data}{A data frame which defines the coordinates and the environmental variables for which the predicted intensity will be calculated. By default, this uses the coordinates and environmental variables from the object supplied to the \code{fit} argument.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates,
as in \code{\link{sampleQuad}}.}
  \item{asp}{Aspect of the plot, with "iso" as default. See the documentation for the \code{levelplot} function of the \code{lattice} pacakge for more details.}
  \item{ylab}{Label for the \eqn{y}-axis. Blank by default.}
  \item{xlab}{Label for the \eqn{x}-axis. Blank by default}
  \item{col.regions}{A vector of colours to define the intensity gradient. See the documentation for the \code{levelplot} function of the \code{lattice} pacakge for more details.}
  \item{cuts}{The number of levels the color gradient for the intensity surface is divided into. See the documentation for the \code{levelplot} function of the \code{lattice} pacakge for more details.}
  \item{cex}{Character size for axis labels.}
  \item{main.text}{Title of the plot.} 
  \item{cex.color}{Character size for the colorkey labels.}
}
\details{This function will compute the predicted intensity of a fitted \code{ppmlasso} object using the model within the regularisation path which optimises the criterion specified in the call to \code{\link{ppmlasso}}.}
\author{Ian W. Renner}
\examples{data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
ppm.form = ~poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2) + poly(D_MAIN_RDS, D_URBAN, degree = 2)
ppm.fit  = ppmlasso(ppm.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, n.fits = 20,
writefile = FALSE)
plotFit(ppm.fit)
}

