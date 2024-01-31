\name{makeMask}
\alias{makeMask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeMask(dat.ppm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat.ppm}{
%%     ~~Describe \code{dat.ppm} here~~
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
function (dat.ppm) 
{
    if (is(dat.ppm, "list")) {
        dat.ppm = dat.ppm$X
    }
    dat.q = dat.ppm[dat.ppm$Pres == 0, ]
    ux = sort(unique(dat.q$X))
    uy = sort(unique(dat.q$Y))
    nx = length(ux)
    ny = length(uy)
    col.ref = match(dat.q$X, ux)
    row.ref = match(dat.q$Y, uy)
    all.vec = rep(0, max(row.ref) * max(col.ref))
    vec.ref = (col.ref - 1) * max(row.ref) + row.ref
    all.vec[vec.ref] = 1
    mask.out = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, 
        ux))
    mask.out
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
