\name{ppmlasso-internal}
\alias{blocks}
\alias{calculateGCV}
\alias{calculateLikelihood}
\alias{calculateUnpenalisedLikelihood}
\alias{catConvert}
\alias{catFrame}
\alias{decimalCount}
\alias{etaFromMu}
\alias{interp}
\alias{irlsUpdate}
\alias{makeMask}
\alias{muFromEta}
\alias{polyNames}
\alias{ppmSpatstat}
\alias{scoreIntercept}
\alias{singleLasso}
\alias{spatRes}
\alias{standardiseX}
\alias{weights}
\alias{zapCoord}
\title{Internal ppmlasso functions}
\description{These are internal functions called by the main functions of \code{ppmlasso}, but are not generally used directly by the user.}
\usage{
blocks(n.blocks, block.scale, dat, seed = 1)
calculateGCV(y, X, ob.wt, b.lasso, lambda, alpha = alpha, 
unp.likelihood = unp.likelihood, penalty = TRUE, family = "poisson", 
mu.min = 1.e-16, mu.max = 1.e16, eta.min = log(1.e-16), eta.max = log(1.e16),
tol = 1.e-9, area.int = FALSE)
calculateLikelihood(X, family, ob.wt, mu, y, alpha, lambda, beta, penalty = FALSE)
calculateUnpenalisedLikelihood(ob.wt, y, mu)
catConvert(env.frame)
catFrame(cat.mat)
decimalCount(x, max.dec = max(10, max(nchar(x))), tol = 0.1)
etaFromMu(mu, family, mu.min = 1.e-16, mu.max = 1/mu.min, eta.min, eta.max)
interp(sp.xy, sp.scale, f, back.xy, coord = c("X", "Y"))
irlsUpdate(y, X, ob.wt, is.in, signs, eta, mu, alpha, lambda, beta.old, 
penalty = FALSE, family, mu.min = 1.e-16, mu.max = 1/mu.min, eta.min, 
eta.max, tol = tol)
makeMask(dat.ppm)
muFromEta(eta, family, mu.min = 1.e-16, mu.max = 1/mu.min)
polyNames(X)
ppmSpatstat(fit)
scoreIntercept(y, X.des, ob.wt = rep(1, length(y)), area.int = FALSE, int = NA, family)
singleLasso(y, X, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA,
intercept = NA, family = "gaussian", tol = 1.e-9, gamma = 0, init.coef = NA, 
mu.min = 1.e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25,
standardise = TRUE)
spatRes(env.grid, coord = c("X", "Y"))
standardiseX(mat)
weights(sp.xy, quad.xy, coord = c("X", "Y"))
zapCoord(x, numdec = decimalCount(x))
}
\author{Ian W. Renner}
\details{These are generally not called by the user.}
\keyword{internal}