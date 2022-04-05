\name{envelope.ppmlasso}
\alias{envelope.ppmlasso}
\title{Calculates simulation envelopes for goodness-of-fit}
\description{This function is analogous to the \code{envelope} function of the \code{spatstat} package.}
\usage{
\method{envelope}{ppmlasso}(Y, fun = Kest, \dots)
}
\arguments{
  \item{Y}{A fitted regularisation path of point process models. The simulation envelopes will be calculated for the model that optimises the given criterion.}
  \item{fun}{The summary function to be computed for the given point process model. See the help file for the \code{envelope} function of the \code{spatstat} package for more details.}
  \item{\dots}{Other arguments for producing diagnostic plots, as given by the \code{envelope} function of the \code{spatstat} package.}
}
\details{See the help file for \code{envelope} in the \code{spatstat} package for further details of simulation envelopes.}
\references{Baddeley, A.J. & Turner, R. (2005). Spatstat: an R package for analyzing spatial
point patterns. \emph{Journal of Statistical Software} \bold{12}, 1-42.}
\author{Ian W. Renner}
\seealso{
\code{\link{diagnose.ppmlasso}}, for residual plots inherited from \code{spatstat}.
}
\examples{
data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
ppm.form = ~poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2) + poly(D_MAIN_RDS, D_URBAN, degree = 2)
ppm.fit  = ppmlasso(ppm.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, n.fits = 20,
writefile = FALSE)
envelope(ppm.fit, Kinhom, nsim = 20)
}