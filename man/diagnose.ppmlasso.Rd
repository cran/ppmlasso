\name{diagnose.ppmlasso}
\alias{diagnose.ppmlasso}
\alias{diagnose}
\title{Create diagnostic plots for a fitted point process model.}
\description{This function is analogous to the \code{diagnose.ppm} function of the \code{spatstat} package.}
\usage{
\method{diagnose}{ppmlasso}(object, \dots)
}
\arguments{
  \item{object}{A fitted regularisation path of point process models. The diagnostic plots will be created for the model that optimises the given criterion.}
  \item{\dots}{Other arguments for producing diagnostic plots, as given by the \code{diagnose.ppm} function of the \code{spatstat} package.}
}
\details{See the help file for \code{diagnose.ppm} in the \code{spatstat} package for further details of diagnostic plots.}
\references{Baddeley, A.J. & Turner, R. (2005). Spatstat: an R package for analyzing spatial
point patterns. \emph{Journal of Statistical Software} \bold{12}, 1-42.}
\author{Ian W. Renner}
\seealso{
\code{\link{envelope.ppmlasso}}, for other goodness-of-fit functions inherited from \code{spatstat}.
}
\examples{
data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
ppm.form = ~poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2) + poly(D_MAIN_RDS, D_URBAN, degree = 2)
ppm.fit  = ppmlasso(ppm.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, n.fits = 20,
writefile = FALSE)
diagnose(ppm.fit, which = "smooth", type = "Pearson")
}
\keyword{residuals}
