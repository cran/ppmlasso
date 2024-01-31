\name{spatRes}
\alias{spatRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spatRes(env.grid, coord = c("X", "Y"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env.grid}{
%%     ~~Describe \code{env.grid} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
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
function (env.grid, coord = c("X", "Y")) 
{
    x.col = which(names(env.grid) == coord[1])
    y.col = which(names(env.grid) == coord[2])
    x.uq = sort(unique(env.grid[, x.col]))
    y.uq = sort(unique(env.grid[, y.col]))
    n.dec = max(unlist(lapply(x.uq, decimalCount)))
    x.diff = diff(x.uq)
    y.diff = diff(y.uq)
    x.dec = unlist(lapply(x.diff, decimalCount))
    y.dec = unlist(lapply(y.diff, decimalCount))
    x.step = min(floor(x.diff * 10^max(x.dec) + 0.1))/(10^max(x.dec))
    y.step = min(floor(y.diff * 10^max(y.dec) + 0.1))/(10^max(y.dec))
    return(list(x.step = x.step, y.step = y.step))
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
