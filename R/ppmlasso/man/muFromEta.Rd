\name{muFromEta}
\alias{muFromEta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
muFromEta(eta, family, mu.min = 1e-16, mu.max = 1/mu.min)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{
%%     ~~Describe \code{eta} here~~
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
function (eta, family, mu.min = 1e-16, mu.max = 1/mu.min) 
{
    mu = family$linkinv(eta)
    mu[mu < mu.min] = mu.min
    mu[mu > mu.max] = mu.max
    mu
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
