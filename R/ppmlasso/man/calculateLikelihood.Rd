\name{calculateLikelihood}
\alias{calculateLikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateLikelihood(X, family, ob.wt, mu, y, alpha, lambda, beta, penalty = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{ob.wt}{
%%     ~~Describe \code{ob.wt} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{penalty}{
%%     ~~Describe \code{penalty} here~~
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
function (X, family, ob.wt, mu, y, alpha, lambda, beta, penalty = FALSE) 
{
    if (penalty == FALSE) {
        lambda = rep(0, dim(X)[2])
    }
    if (family$family == "poisson") {
        like = sum(ob.wt * (y * log(mu) - mu)) - sum(log(1:sum(y > 
            0))) - sum(alpha * as.vector(lambda) * abs(beta)) - 
            sum(0.5 * (1 - alpha) * as.vector(lambda) * beta^2)
    }
    if (family$family == "binomial") {
        like = sum(ob.wt * (y * log(mu) + (1 - y) * log(1 - mu))) - 
            sum(as.vector(lambda) * abs(beta))
    }
    if (family$family == "gaussian") {
        like = sum(ob.wt * (y * mu - 0.5 * mu^2 - 0.5 * y^2 - 
            0.5 * log(2 * pi))) - sum(as.vector(lambda) * abs(beta))
    }
    like
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
