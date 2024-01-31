\name{calculateGCV}
\alias{calculateGCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateGCV(y, X, ob.wt, b.lasso, lambda, alpha = alpha, unp.likelihood = unp.likelihood, penalty = TRUE, family = "poisson", mu.min = 1e-16, mu.max = 1e+16, eta.min = log(1e-16), eta.max = log(1e+16), tol = 1e-09, area.int = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{ob.wt}{
%%     ~~Describe \code{ob.wt} here~~
}
  \item{b.lasso}{
%%     ~~Describe \code{b.lasso} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{unp.likelihood}{
%%     ~~Describe \code{unp.likelihood} here~~
}
  \item{penalty}{
%%     ~~Describe \code{penalty} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{mu.min}{
%%     ~~Describe \code{mu.min} here~~
}
  \item{mu.max}{
%%     ~~Describe \code{mu.max} here~~
}
  \item{eta.min}{
%%     ~~Describe \code{eta.min} here~~
}
  \item{eta.max}{
%%     ~~Describe \code{eta.max} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{area.int}{
%%     ~~Describe \code{area.int} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, X, ob.wt, b.lasso, lambda, alpha = alpha, unp.likelihood = unp.likelihood, 
    penalty = TRUE, family = "poisson", mu.min = 1e-16, mu.max = 1e+16, 
    eta.min = log(1e-16), eta.max = log(1e+16), tol = 1e-09, 
    area.int = FALSE) 
{
    is.in = abs(b.lasso) > 1e-07
    signs = sign(b.lasso)
    eta = X \%*\% b.lasso
    mu = exp(eta)
    n.pres = sum(y > 0)
    wt = irlsUpdate(y, X, ob.wt, is.in, signs, eta, mu, alpha, 
        lambda, b.lasso, penalty = TRUE, family = family, mu.min, 
        mu.max, eta.min, eta.max, tol)$wt
    xwx = wt \%*\% X[, is.in]
    bpinv = rep(0, length(b.lasso))
    bpinv[is.in] = 1/abs(b.lasso[is.in])
    ginv = diag(bpinv)
    if (sum(is.in) > 1) {
        eff.df = sum(diag(solve(xwx + diag(as.vector(lambda[is.in])) \%*\% 
            ginv[is.in, is.in]) \%*\% xwx)) + area.int
    }
    if (sum(is.in) == 1) {
        eff.df = sum(diag(solve(xwx + diag(as.matrix(lambda[is.in])) \%*\% 
            ginv[is.in, is.in]) \%*\% xwx)) + area.int
    }
    if (family$family == "poisson") {
        dev = -2 * (sum(ob.wt * (y * log(mu) - mu)) + sum(y > 
            0))
    }
    if (family$family == "binomial") {
        dev = -2 * unp.likelihood
    }
    if (family$family == "gaussian") {
        dev = -2 * (unp.likelihood - calculateLikelihood(X, family, 
            ob.wt[y > 0], y[y > 0], y, alpha, lambda, beta = b.lasso, 
            penalty = FALSE))
    }
    if (eff.df < n.pres) {
        gcv = dev/(n.pres * (1 - (eff.df)/n.pres)^2)
    }
    if (eff.df >= n.pres) {
        gcv = NA
    }
    return(list(gcv = gcv, dev = dev, eff.df = eff.df))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
