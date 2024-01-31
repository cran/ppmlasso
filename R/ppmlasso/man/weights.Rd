\name{weights}
\alias{weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
weights(sp.xy, quad.xy, coord = c("X", "Y"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.xy}{
%%     ~~Describe \code{sp.xy} here~~
}
  \item{quad.xy}{
%%     ~~Describe \code{quad.xy} here~~
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
function (sp.xy, quad.xy, coord = c("X", "Y")) 
{
    sp.col = c(which(names(sp.xy) == coord[1]), which(names(sp.xy) == 
        coord[2]))
    quad.col = c(which(names(quad.xy) == coord[1]), which(names(quad.xy) == 
        coord[2]))
    X.inc = sort(unique(quad.xy[, quad.col[1]]))[2] - sort(unique(quad.xy[, 
        quad.col[1]]))[1]
    Y.inc = sort(unique(quad.xy[, quad.col[2]]))[2] - sort(unique(quad.xy[, 
        quad.col[2]]))[1]
    quad.0X = min(quad.xy[, quad.col[1]]) - floor(min(quad.xy[, 
        quad.col[1]])/X.inc) * X.inc
    quad.0Y = min(quad.xy[, quad.col[2]]) - floor(min(quad.xy[, 
        quad.col[2]])/Y.inc) * Y.inc
    X = c(sp.xy[, quad.col[1]], quad.xy[, quad.col[1]])
    Y = c(sp.xy[, quad.col[2]], quad.xy[, quad.col[2]])
    round.X = round((X - quad.0X)/X.inc) * X.inc
    round.Y = round((Y - quad.0Y)/Y.inc) * Y.inc
    round.id = paste(round.X, round.Y)
    round.table = table(round.id)
    wt = X.inc * Y.inc/as.numeric(round.table[match(round.id, 
        names(round.table))])
    wt
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
