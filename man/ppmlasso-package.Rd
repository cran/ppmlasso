\name{ppmlasso-package}
\alias{ppmlasso-package}
\docType{package}
\title{PPM-LASSO: Point process models with LASSO-type penalties}
\description{This package contains tools to fit point process models with sequences of LASSO penalties ("regularisation paths"). 
Regularisation paths of Poisson point process models or area-interaction models can be fitted with LASSO, 
adaptive LASSO or elastic net penalties. A number of criteria are available to judge the bias-variance tradeoff.}
\details{
The key functions in ppmlasso are as follows:

\bold{Useful pre-analysis functions:}
	\describe{
		\item{\code{\link{findRes}}}{Determine the optimal spatial resolution at which to perform analysis}
		\item{\code{\link{getEnvVar}}}{Interpolate environmental data to species presence locations}
		\item{\code{\link{griddify}}}{Ensure a matrix of environmental data is on a rectangular grid}
		\item{\code{\link{ppmdat}}}{Calculate observation weights and set up design matrix for fitting}
		\item{\code{\link{pointInteractions}}}{Calculate interpoint interactions for fitting area-interaction models}
		\item{\code{\link{sampleQuad}}}{Set up a regular grid of quadrature points}
	}
\bold{Creating regularisation paths of point process models:}
	\describe{
		\item{\code{\link{ppmlasso}}}{Fit a regularisation path of point process models}
		\item{\code{\link{plotFit}}}{Plot the fitted intensity of a \code{ppmlasso} object}
		\item{\code{\link{plotPath}}}{Plot the regularisation path of a \code{ppmlasso} object}
		\item{\code{\link{print.ppmlasso}}}{Print output from a \code{ppmlasso} object}
		\item{\code{\link{predict.ppmlasso}}}{Make predictions from a fitted point process model to new data}
	}
\bold{Checking assumptions:}
	\describe{
		\item{\code{\link{diagnose.ppmlasso}}}{Create diagnostic residual plots of \code{ppmlasso} object}
		\item{\code{\link{envelope.ppmlasso}}}{Create simulation envelope for goodness-of-fit checks on a \code{ppmlasso} object}
	}
}
\author{Ian W. Renner

Maintainer: Ian W. Renner <Ian.Renner@newcastle.edu.au>
}
\references{
Renner, I.W. & Warton, D.I. (2013). Equivalence of MAXENT and Poisson point process models for 
species distribution modeling in ecology \emph{Biometrics} \bold{69}, 274-281.

Renner, I.W. et al (2015). Point process models for presence-only analysis. \emph{Methods in Ecology and Evolution} \bold{6}, 366-379.

Renner, I.W., Warton, D.I., & Hui, F.K.C. (2021). What is the effective sample size of a spatial point process? \emph{Australian & New Zealand Journal of Statistics} \bold{63}, 144-158.

}
\keyword{ package }
\examples{
# Fit a regularisation path of Poisson point process models
data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
ppm.form = ~ poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2, raw = TRUE)
ppm.fit  = ppmlasso(ppm.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, n.fits = 20,
writefile = FALSE)

# Fit a regularisation path of area-interaction models
data(BlueMountains)
ai.form  = ~ poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2, raw = TRUE)
ai.fit   = ppmlasso(ai.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, 
family = "area.inter", r = 2, availability = BlueMountains$availability, n.fits = 20,
writefile = FALSE)

# Print a ppmlasso object
print(ppm.fit, out = "model")

# Residual plot of a ppmlasso object
diagnose(ppm.fit, which = "smooth", type = "Pearson")

# Make predictions
pred.mu = predict(ppm.fit, newdata = sub.env)

# Plot the intensity from a fitted ppmlasso object
plotFit(ppm.fit)

# Plot the regularisation path from a fitted ppmlasso object
plotPath(ppm.fit)

}
