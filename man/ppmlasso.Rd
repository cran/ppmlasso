\name{ppmlasso}
\alias{ppmlasso}
\title{Fit point process models with LASSO penalties}
\description{The ppmlasso function fits point process models (either Poisson or area-interaction
models) with a sequence of LASSO, adaptive LASSO or elastic net penalties (a "regularisation path").}
\usage{
ppmlasso(formula, sp.xy, env.grid, sp.scale, coord = c("X", "Y"), 
data = ppmdat(sp.xy = sp.xy, sp.scale = sp.scale, back.xy = env.grid, coord = c("X","Y"),
sp.file = NA, quad.file = NA, datfilename = "PPMDat", writefile = writefile), lamb = NA,
n.fits = 200, ob.wt = NA, criterion = "bic", alpha = 1, family = "poisson", tol = 1.e-9, 
gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, r = NA, interactions = NA, 
availability = NA, max.it = 25, min.lamb = -10, standardise = TRUE, n.blocks = NA, 
block.size = sp.scale*100, seed = 1, writefile = TRUE)
}
\arguments{
  \item{formula}{The formula of the fitted model. For a point process model, the correct form is \code{~ variables}.}
  \item{sp.xy}{A matrix of species locations containing at least one column representing
longitude and one column representing latitude. Environmental variables are
interpolated to the locations of \code{sp.xy} using the \code{\link{getEnvVar}} function, unless the
\code{data} argument is supplied.}
  \item{env.grid}{The geo-referenced matrix of environmental grids. This matrix is used to
generate quadrature points using the \code{\link{sampleQuad}} function, interpolate environmental
data to the species locations of \code{sp.xy} using the \code{\link{getEnvVar}} function,
and calculate observation weights using the \code{\link{ppmdat}} function, unless the \code{data}
argument is supplied. This creates a data matrix \code{data} which provides the
variables for the \code{formula} argument.}
  \item{sp.scale}{The spatial resolution at which to define the regular grid of quadrature
points. \code{\link{sampleQuad}} will subsample from the rows of \code{data} that coincide with
a regular grid at a resolution of \code{sp.scale}.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates.}
  \item{data}{An optional data matrix generated from the \code{\link{ppmdat}} function. Supplying a
matrix to \code{data} is an alternative way of providing the environmental variables
used in the \code{formula} argument, instead of specifying \code{sp.xy} and \code{env.grid}.}
  \item{lamb}{A vector of penalty values that will be used to create the regularisation path.
If \code{lamb = NA}, the penalty values are automatically generated from the
\code{data} and the \code{n.fits} argument.}
  \item{n.fits}{The number of models fitted in the regularisation path. If \code{lamb = NA}, the
\code{n.fits} penalty values will be equally spaced on a logarithmic scale from \eqn{e^{-10}}{exp(-10)}
to \eqn{\lambda_{max}}{\lambda_{max}}, the smallest penalty that shrinks all parameter coefficients to zero.}
  \item{ob.wt}{Quadrature weights, usually inherited from the \code{\link{ppmdat}} function.}
  \item{criterion}{The penalisation criteria to be optimised by the regularisation path. The
options include \code{"aic"}, \code{"bic"}, \code{"blockCV"}, \code{"hqc"}, \code{"gcv"}, \code{"nlgcv"} and \code{"msi"}.}
  \item{alpha}{The elastic net parameter. The form of the penalty is \deqn{\alpha*\lambda*\sum_{j = 1}^p |\beta_j| + (1 - \alpha)*\lambda*\sum_{j = 1}^p (\beta_j)^2.}{\alpha*\lambda*\sum_{j = 1}^p |\beta_j| + (1 - \alpha)*\lambda*\sum_{j = 1}^p (\beta_j)^2.} The default value \code{alpha = 1} corresponds to a LASSO penalty,
while \code{alpha = 0} corresponds to a ridge regression penalty.}
  \item{family}{The family of models to be fitted -- \code{family = "poisson"} for Poisson point process models
or \code{family = "area.inter"} for area-interaction models.}
  \item{tol}{The convergence threshold for the descent algorithm. The algorithm continues
for a maximum of \code{max.it} iterations until the difference in likelihood between
successive fits falls below \code{tol}.}
  \item{gamma}{The exponent of the adaptive weights for the adaptive LASSO penalty. The
default value \code{gamma = 0} corresponds to a normal LASSO penalty.}
  \item{init.coef}{The initial coefficients used for an adaptive LASSO penalty.}
  \item{mu.min}{The threshold for small fitted values. Any fitted value less than the threshold
is set to \code{mu.min}.}
  \item{mu.max}{The threshold for large fitted values. Any fitted value larger than the threshold
will be set to \code{mu.max}.}
  \item{r}{The radius of point interactions, required if \code{family = "area.inter"}.}
  \item{interactions}{A vector of point interactions calculated from the \code{\link{pointInteractions}}
function necessary for fitting area-interaction models. If \code{interactions = NA}
and \code{family = "area.inter"}, point interactions will be automatically calculated
for radius \code{r} to the locations of \code{data}.}
  \item{availability}{An optional binary matrix used in calculating point interactions indicating
whether locations are available (1) or not (0). See \code{\link{pointInteractions}} for more details.}
  \item{max.it}{The maximum number of iterations of the descent algorithm for fitting the
model.}
  \item{min.lamb}{The power \eqn{x} of smallest penalty \eqn{e^x} among the \code{n.fits} models.}
  \item{standardise}{A logical argument indicating whether the environmental variables
should be standardised to have mean 0 and variance 1. It is recommended that variables are
standardised for analysis.}
  \item{n.blocks}{This argument controls the number of cross validation groups into which the spatial blocks
are divided if the \code{criterion} argument is set to \code{"blockCV"}. See details.}
  \item{block.size}{The length of the edges for the spatial blocks created if the \code{criterion} argument 
is set to \code{"blockCV"}. Only square spatial blocks are currently supported. See details.}
  \item{seed}{The random seed used for controlling the allocation of spatial blocks to cross validation groups
if the \code{criterion} argument is set to \code{"blockCV"}.}
  \item{writefile}{A logical argument passed to the \code{\link{ppmdat}} function to determine whether its output should be written to a file or not, set to \code{TRUE} by default. See the documentation for \code{\link{ppmdat}} for details.}
}
\details{This function fits a regularisation path of point process models provided a list of species locations
and a geo-referenced grid of environmental data. It is assumed that Poisson point process models (Warton & Shepherd, 2010)
fit intensity as a log-linear model of environmental covariates, and that area-interaction models
(Widom & Rowlinson, 1970; Baddeley & van Lieshout, 1995) fit conditional intensity as a log-linear model of
environmental covariates and point interactions. Parameter coefficients are estimated by maximum likelihood for
Poisson point process models and by maximum pseudolikelihood (Besag, 1977) for area-interaction models. The expressions
for both the likelihood and pseudolikelihood involve an intractable integral which is approximated using a quadrature
scheme (Berman & Turner, 1992).

Each model in the regularisation path is fitted by extending the Osborne descent algorithm (Osborne, 2000) to
generalised linear models with penalised iteratively reweighted least squares.

Three classes of penalty \eqn{p(\beta)} are available for the vector of parameter coefficients \eqn{\beta}:

For the LASSO (Tibshirani, 1996), \eqn{p(\beta) = \lambda*\sum_{j = 1}^p |\beta_j|}

For the adaptive LASSO (Zou, 2006), \eqn{p(\beta) = \lambda*\sum_{j = 1}^p w_j*|\beta_j|}, where \eqn{w_j = 1/|\hat{\beta}_{init, j}|^\gamma}
for some initial estimate of parameters \eqn{\hat{\beta}_{init}}.

For the elastic net (Zou & Hastie, 2005), \eqn{\alpha*\lambda*\sum_{j = 1}^p |\beta_j| + (1 - \alpha)*\lambda*\sum_{j = 1}^p (\beta_j)^2}{\alpha*\lambda*\sum_{j = 1}^p |\beta_j| + (1 - \alpha)*\lambda*\sum_{j = 1}^p (\beta_j)^2}.
Note that this form of the penalty is a restricted case of the general elastic net penalty.

There are various criteria available for managing the bias-variance tradeoff (Renner, 2013).
The default choice is BIC, the Bayesian Information Criterion, which has been shown to have
good performance. 

An alternative criterion useful when data are sparse is MSI, the maximum score 
of the intercept model (Renner, in prep). For a set of \eqn{m} presence locations, the
MSI penalty is \eqn{\lambda_{MSI} = \lambda_{max}/\sqrt{m}}, where \eqn{\lambda_{max}} is the smallest penalty
that shrinks all environmental coefficients to zero. The MSI penalty differs from the other criteria in that does
not require an entire regularisation path to be fitted.

It is also possible to control the magnitude of the penalty by spatial cross validation by setting the 
\code{criterion} argument to \code{"blockCV"}. The study region is then divided into square blocks with edge
lengths controlled by the \code{block.size} argument, which are assigned to one of a number of cross validation
groups controlled by the \code{n.groups} argument. The penalty which maximises the predicted log-likelihood is
chosen.
}
\value{An object of class \code{"ppmlasso"}, with elements:
\item{betas}{A matrix of fitted coefficients of the \code{n.fits} models.}
\item{lambdas}{A vector containing the \code{n.fits} penalty values.}
\item{likelihoods}{A vector containing the likelihood of \code{n.fits} fitted models.}
\item{pen.likelihoods}{A vector containing the penalised likelihood of \code{n.fits} fitted
models.}
\item{beta}{A vector containing the coefficients of the model that optimises the specified \code{criterion}.}
\item{lambda}{The penalty value of the model that optimises the specified \code{criterion}.}
\item{mu}{A vector of fitted values from the model that optimises the specified \code{criterion}.}
\item{likelihood}{The likelihood of the model that optimises the specified \code{criterion}.}
\item{criterion}{The specified \code{criterion} of the function call.}
\item{family}{The specified \code{family} of the function call.}
\item{gamma}{The specified \code{gamma} of the function call.}
\item{alpha}{The specified \code{alpha} of the function call.}
\item{init.coef}{The specified \code{init.coef} of the function call.}
\item{criterion.matrix}{A matrix with \code{n.fits} rows corresponding to the observed values
of AIC, BIC, HQC, GCV, and non-linear GCV.}
\item{data}{The design matrix. For the point process models fitted with this function,
\code{mu = e^{data*beta}}.}
\item{pt.interactions}{The calculated point interactions.}
\item{wt}{The vector of quadrature weights.}
\item{pres}{A vector indicating presence (1) or quadrature point (0).}
\item{x}{A vector of point longitudes.}
\item{y}{A vector of point latitudes.}
\item{r}{The radius of point interactions.}
\item{call}{The function call.}
\item{formula}{The \code{formula} argument.}
\item{s.means}{If \code{standardise = TRUE}, the means of each column of \code{data} prior to standardisation.}
\item{s.sds}{If \code{standardise = TRUE}, the standard deviations of each column of \code{data} prior to standardisation.}
\item{cv.group}{The cross validation group associated with each point in the data set.}
\item{n.blocks}{The number of cross validation groups specified.}
}
\references{
Baddeley, A.J. & van Lieshout, M.N.M. (1995). Area-interaction point processes.
\emph{Annals of the Institute of Statistical Mathematics} \bold{47}, 601-619.

Berman, M. & Turner, T.R. (1992). Approximating point process likelihoods with
GLIM. \emph{Journal of the Royal Statistics Society, Series C} \bold{41}, 31-38.

Besag, J. (1977). Some methods of statistical analysis for spatial data. \emph{Bulletin of
the International Statistical Institute} \bold{47}, 77-91.

Osborne, M.R., Presnell, B., & Turlach, B.A. (2000). On the lasso and its dual.
\emph{Journal of Computational and Graphical Statistics} \bold{9}, 319-337.

Renner, I.W. & Warton, D.I. (2013). Equivalence of MAXENT and Poisson point process 
models for species distribution modeling in ecology. \emph{Biometrics} \bold{69}, 274-281.

Renner, I.W. (2013). Advances in presence-only methods in ecology.
\url{https://unsworks.unsw.edu.au/fapi/datastream/unsworks:11510/SOURCE01}

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. \emph{Journal of
the Royal Statistical Society, Series B} \bold{58}, 267-288.

Warton, D.I. & Shepherd, L.C. (2010). Poisson point process models solve the
"pseudo-absence problem" for presence-only data in ecology. \emph{Annals of Applied
Statistics} \bold{4}, 1383-1402.

Widom, B. & Rowlinson, J.S. (1970). New model for the study of liquid-vapor
phase transitions. \emph{The Journal of Chemical Physics} \bold{52}, 1670-1684.

Zou, H. (2006). The adaptive lasso and its oracle properties. \emph{Journal of the American
Statistical Association} \bold{101}, 1418-1429.

Zou, H. & Hastie, T. (2005). Regularization and variable selection via the elastic
net. \emph{Journal of the Royal Statistical Society, Series B} \bold{67}, 301-320.}
\author{Ian W. Renner}
\seealso{
\code{\link{print.ppmlasso}} for printing features of the fitted regularisation path.

\code{\link{predict.ppmlasso}} for predicting intensity for a set of new data.

\code{\link{envelope.ppmlasso}} for constructing a K-envelope of the model which optimises
the given criterion from the \code{spatstat} package.

\code{\link{diagnose.ppmlasso}} for diagnostic plots from the \code{spatstat} package.
}
\examples{
# Fit a regularisation path of Poisson point process models
data(BlueMountains)
sub.env = BlueMountains$env[BlueMountains$env$Y > 6270 & BlueMountains$env$X > 300,]
sub.euc = BlueMountains$eucalypt[BlueMountains$eucalypt$Y > 6270 & BlueMountains$eucalypt$X > 300,]
ppm.form = ~ poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2)
ppm.fit  = ppmlasso(ppm.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, n.fits = 20,
writefile = FALSE)

#Fit a regularisation path of area-interaction models
data(BlueMountains)
ai.form  = ~ poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2)
ai.fit   = ppmlasso(ai.form, sp.xy = sub.euc, env.grid = sub.env, sp.scale = 1, 
family = "area.inter", r = 2, availability = BlueMountains$availability, n.fits = 20,
writefile = FALSE)
}
\keyword{LASSO}
\keyword{Point process models}
\keyword{Regularisation path}
