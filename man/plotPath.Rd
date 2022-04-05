\name{plotPath}
\alias{plotPath}
\title{Plot of the regularisation path of a \code{ppmlasso} model}
\description{This function produces a trace plot of the coefficient estimates of a fitted \code{ppmlasso} object and identifies the models which optimise the various penalisation criteria.}
\usage{
plotPath(fit, colors = c("gold", "green3", "blue", "brown", "pink"), logX = TRUE)
}
\arguments{
  \item{fit}{A fitted \code{ppmlasso} object.}
  \item{colors}{A vector of colours for each criterion: \code{"aic"}, \code{"bic"}, \code{"hqc"}, \code{"gcv"}, and \code{"nlgcv"}.}
  \item{logX}{A logical argument to indicate whether the plot should utilise a logarithmic scale on the x-axis (the default) or not.}
}
\details{A fitted \code{ppmlasso} object contains a matrix called \code{betas} which stores the coefficient estimates of each of the \code{n.fits} models fitted. This function produces a traceplot of these coefficient estimates for each environmental variable and highlights the models which optimise each of the penalisation criteria.}
\author{Ian W. Renner}
\examples{data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
ppm.form = ~poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2) + poly(D_MAIN_RDS, D_URBAN, degree = 2)
ppm.fit  = ppmlasso(ppm.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, n.fits = 20,
writefile = FALSE)
plotPath(ppm.fit)
}