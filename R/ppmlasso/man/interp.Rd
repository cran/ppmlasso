\name{interp}
\alias{interp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interp(sp.xy, sp.scale, f, back.xy, coord = c("X", "Y"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.xy}{
%%     ~~Describe \code{sp.xy} here~~
}
  \item{sp.scale}{
%%     ~~Describe \code{sp.scale} here~~
}
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{back.xy}{
%%     ~~Describe \code{back.xy} here~~
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
function (sp.xy, sp.scale, f, back.xy, coord = c("X", "Y")) 
{
    options(scipen = 999)
    x.dat = sp.xy[, which(names(sp.xy) == coord[1])]
    y.dat = sp.xy[, which(names(sp.xy) == coord[2])]
    x.back = back.xy[, which(names(back.xy) == coord[1])]
    y.back = back.xy[, which(names(back.xy) == coord[2])]
    res = spatRes(back.xy, coord = coord)
    grid = data.table(x.back, y.back, f, key = c("x.back", "y.back"))
    ux = sort(unique(x.back))
    uy = sort(unique(y.back))
    x.col = which(names(back.xy) == coord[1])
    y.col = which(names(back.xy) == coord[2])
    x.step = res$x.step
    y.step = res$y.step
    x.o = min(back.xy[, x.col]) - floor(min(back.xy[, x.col])/x.step) * 
        x.step
    y.o = min(back.xy[, y.col]) - floor(min(back.xy[, y.col])/y.step) * 
        y.step
    x.1 = floor((x.dat - x.o)/sp.scale) * sp.scale + x.o
    y.1 = floor((y.dat - y.o)/sp.scale) * sp.scale + y.o
    x.2 = pmin(x.1 + sp.scale, max(ux))
    y.2 = pmin(y.1 + sp.scale, max(uy))
    w11 = (x.2 - x.dat) * (y.2 - y.dat)/((x.2 - x.1) * (y.2 - 
        y.1))
    w12 = (x.2 - x.dat) * (y.dat - y.1)/((x.2 - x.1) * (y.2 - 
        y.1))
    w21 = (x.dat - x.1) * (y.2 - y.dat)/((x.2 - x.1) * (y.2 - 
        y.1))
    w22 = (x.dat - x.1) * (y.dat - y.1)/((x.2 - x.1) * (y.2 - 
        y.1))
    f11 = grid[list(x.1, y.1)]$f
    f12 = grid[list(x.1, y.2)]$f
    f21 = grid[list(x.2, y.1)]$f
    f22 = grid[list(x.2, y.2)]$f
    c11 = 1 - is.na(f11)
    c12 = 1 - is.na(f12)
    c21 = 1 - is.na(f21)
    c22 = 1 - is.na(f22)
    env.wt.mat = cbind(f11 * w11 * c11, f12 * w12 * c12, f21 * 
        w21 * c21, f22 * w22 * c22)
    f.interp = apply(env.wt.mat, 1, sum, na.rm = TRUE)/(w11 * 
        c11 + w12 * c12 + w21 * c21 + w22 * c22)
    f.interp
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
