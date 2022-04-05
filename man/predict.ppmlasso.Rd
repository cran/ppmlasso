\name{predict.ppmlasso}
\alias{predict.ppmlasso}
\title{Prediction to new data from a fitted regularisation path}
\description{Given a fitted regularisation path produced by \code{\link{ppmlasso}}, this function will predict the intensity for a new
set of data.}
\usage{
\method{predict}{ppmlasso}(object, \dots, newdata, interactions = NA)
}
\arguments{
  \item{object}{A fitted regularisation path produced by \code{\link{ppmlasso}}.}
  \item{\dots}{Additional arguments impacting the prediction calculations.}
  \item{newdata}{A data frame of new environmental data for which predicted intensities are computed.}
  \item{interactions}{A vector of point interactions for predictions of area-interaction models.}
}
\value{A vector of predicted intensities corresponding to the environmental data provided in the \code{newdata} argument.}
\author{Ian W. Renner}
\seealso{
\code{\link{ppmlasso}} for fitting a regularisation path of point process models.
}
\examples{
data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
ppm.form = ~ poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2, raw = TRUE)
ppm.fit  = ppmlasso(ppm.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, n.fits = 20,
writefile = FALSE)
pred.mu  = predict(ppm.fit, newdata = sub.env)
}
\keyword{predict}
