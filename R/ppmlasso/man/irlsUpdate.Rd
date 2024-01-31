\name{irlsUpdate}
\alias{irlsUpdate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
irlsUpdate(y, X, ob.wt, is.in, signs, eta, mu, alpha, lambda, beta.old, penalty = FALSE, family, mu.min = 1e-16, mu.max = 1/mu.min, eta.min, eta.max, tol = tol)
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
  \item{is.in}{
%%     ~~Describe \code{is.in} here~~
}
  \item{signs}{
%%     ~~Describe \code{signs} here~~
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{beta.old}{
%%     ~~Describe \code{beta.old} here~~
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
function (y, X, ob.wt, is.in, signs, eta, mu, alpha, lambda, 
    beta.old, penalty = FALSE, family, mu.min = 1e-16, mu.max = 1/mu.min, 
    eta.min, eta.max, tol = tol) 
{
    if (penalty == FALSE) {
        is.in = rep(TRUE, dim(X)[2])
        lambda = rep(0, dim(X)[2])
        signs = rep(0, dim(X)[2])
    }
    vari = family$variance(mu)
    deriv = 1/vari
    z = eta + (y - mu) * deriv
    weii = ob.wt * 1/(deriv^2 * vari)
    Xw = t(as.vector(weii) * t(t(X[, is.in])))
    Xw.s = t(as.vector(weii) * t(t(X)))
    epsilon = diag(rep(sqrt(tol), sum(is.in)))
    if (sum(is.in) == 1) {
        epsilon = sqrt(tol)
    }
    xwx = Xw \%*\% X[, is.in] + epsilon
    qx = qr(xwx)
    dim.X = dim(xwx)
    if (qx$rank >= dim.X[1]) {
        if (dim.X[1] == 1) {
            beta.new = solve(xwx + (1 - alpha) * lambda[is.in]) \%*\% 
                (Xw \%*\% z - as.matrix(alpha * lambda[is.in] * 
                  signs[is.in]))
        }
        if (dim.X[1] > 1) {
            beta.new = solve(xwx + diag((1 - alpha) * lambda[is.in])) \%*\% 
                (Xw \%*\% z - as.matrix(alpha * lambda[is.in] * 
                  signs[is.in]))
        }
        eta.new = as.matrix(X[, is.in]) \%*\% as.matrix(beta.new)
        eta.new[eta.new < eta.min] = eta.min
        eta.new[eta.new > eta.max] = eta.max
        mu.new = family$linkinv(eta.new)
        mu.new[mu.new < mu.min] = mu.min
        mu.new[mu.new > mu.max] = mu.max
        return(list(mu = mu.new, beta = beta.new, eta = eta.new, 
            wt = Xw, XwX = xwx, s.wt = Xw.s, deriv = deriv, v = vari, 
            error = "None"))
    }
    if (qx$rank < dim.X) {
        error.flag = "Singular matrix"
        return(list(error = error.flag))
    }
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
