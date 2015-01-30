\name{ppmlasso-internal}
\alias{blocks}
\alias{CatConvert}
\alias{CatFrame}
\alias{eta.from.mu}
\alias{gcv.calc}
\alias{interp}
\alias{irls.update}
\alias{like.calc}
\alias{mu.from.eta}
\alias{ppm.ss}
\alias{score.int}
\alias{single.lasso}
\alias{standardise.X}
\alias{unp.likelihood}
\alias{weights}
\title{Internal ppmlasso functions}
\description{These are internal functions called by the main functions of ppmlasso, but are not generally used directly by the user.}
\usage{
blocks(n.blocks, block.scale, dat, seed = 1)
CatConvert(env.frame)
CatFrame(cat.mat)
eta.from.mu(mu, family, mu.min = 1e-16, mu.max = 1/mu.min, eta.min, eta.max)
gcv.calc(y, X, ob.wt, b.lasso, lambda, alpha = alpha, unp.likelihood = unp.likelihood, 
penalty = TRUE, family = "poisson", mu.min = 1e-16, mu.max = 1e+16, 
eta.min = log(1e-16), eta.max = log(1e+16), tol = 1e-09, area.int = FALSE)
interp(sp.xy, sp.scale, f, back.xy, coord = c("X", "Y"))
irls.update(y, X, ob.wt, is.in, signs, eta, mu, alpha, lambda, beta.old, 
penalty = FALSE, family, mu.min = 1e-16, mu.max = 1/mu.min, eta.min, eta.max, tol = tol)
like.calc(X, family, ob.wt, mu, y, alpha, lambda, beta, penalty = FALSE)
mu.from.eta(eta, family, mu.min = 1e-16, mu.max = 1/mu.min)
ppm.ss(fit)
score.int(y, X.des, ob.wt = rep(1, length(y)), area.int = FALSE, int = NA, family)
single.lasso(y, X, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA, 
intercept = NA, family = "gaussian", tol = 1e-09, gamma = 0, init.coef = NA, 
mu.min = 1e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25, 
standardise = TRUE)
standardise.X(mat)
unp.likelihood(ob.wt, y, mu)
weights(sp.xy, quad.xy, coord = c("X", "Y"))
}
\author{Ian W. Renner}
\details{These are generally not called by the user.}
\keyword{internal}