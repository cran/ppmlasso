\name{scoreIntercept}
\alias{scoreIntercept}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreIntercept(y, X.des, ob.wt = rep(1, length(y)), area.int = FALSE, int = NA, family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{X.des}{
%%     ~~Describe \code{X.des} here~~
}
  \item{ob.wt}{
%%     ~~Describe \code{ob.wt} here~~
}
  \item{area.int}{
%%     ~~Describe \code{area.int} here~~
}
  \item{int}{
%%     ~~Describe \code{int} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
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
function (y, X.des, ob.wt = rep(1, length(y)), area.int = FALSE, 
    int = NA, family) 
{
    if (area.int == FALSE) {
        int.mod = glm(y ~ 1, family = family, weights = ob.wt)
    }
    if (area.int != FALSE) {
        int.mod = glm(y ~ 1 + int, family = family, weights = ob.wt)
    }
    mu = int.mod$fitted
    vari = family$variance(mu)
    deriv = 1/vari
    wt.mat = ob.wt * 1/(deriv^2 * vari)
    Xw.s = t(as.vector(wt.mat) * t(t(X.des)))
    score = t(as.vector(deriv) * t(Xw.s)) \%*\% (y - mu)
    score
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
