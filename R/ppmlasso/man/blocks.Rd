\name{blocks}
\alias{blocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
blocks(n.blocks, block.scale, dat, seed = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.blocks}{
%%     ~~Describe \code{n.blocks} here~~
}
  \item{block.scale}{
%%     ~~Describe \code{block.scale} here~~
}
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
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
function (n.blocks, block.scale, dat, seed = 1) 
{
    if (is(dat, "list")) {
        dat = dat$X
    }
    cell.group = rep(0, length(dat$X))
    n.groups = ceiling(c(max((dat$X - min(dat$X))/block.scale), 
        max((dat$Y - min(dat$Y))/block.scale)))
    xq = block.scale * (1:n.groups[1]) + min(dat$X)
    for (i.group in 1:n.groups[1]) {
        cell.group = cell.group + as.numeric(dat$X > xq[i.group])
    }
    yq = block.scale * (1:n.groups[2]) + min(dat$Y)
    for (i.group in 1:n.groups[2]) {
        cell.group = cell.group + n.groups[1] * as.numeric(dat$Y > 
            yq[i.group])
    }
    block.group = factor(cell.group)
    set.seed(seed)
    levs = rep(1:n.blocks, length = length(levels(block.group)))
    lev.sample = sample(levs)
    levels(block.group) = lev.sample
    block.group
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
