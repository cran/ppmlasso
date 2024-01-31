\name{griddify}
\alias{griddify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
griddify(envframe, tol = 0.01, coord = c("X", "Y"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{envframe}{
%%     ~~Describe \code{envframe} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
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
function (envframe, tol = 0.01, coord = c("X", "Y")) 
{
    xcol = match(coord[1], names(envframe))
    ycol = match(coord[2], names(envframe))
    quadx = envframe[, xcol]
    quady = envframe[, ycol]
    dx = quadx[duplicated(quadx)]
    dy = quady[duplicated(quady)]
    udx = sort(unique(dx))
    udy = sort(unique(dy))
    diffx = diff(udx)
    diffy = diff(udy)
    x_grain = median(diffx)
    y_grain = median(diffy)
    x_first = udx[which(diffx == x_grain)[1]]
    y_first = udy[which(diffy == y_grain)[1]]
    x_adj = 0 - (x_first - 0.5 * x_grain - round((x_first - 0.5 * 
        x_grain)/x_grain) * x_grain) - 0.5 * x_grain
    y_adj = 0 - (y_first - 0.5 * y_grain - round((y_first - 0.5 * 
        y_grain)/y_grain) * y_grain) - 0.5 * y_grain
    x_fuzzy_high = quadx + x_adj - round((quadx + x_adj)/x_grain) * 
        x_grain + tol * x_grain
    x_fuzzy_low = quadx + x_adj - round((quadx + x_adj)/x_grain) * 
        x_grain - tol * x_grain
    x_in = x_fuzzy_high >= 0 & x_fuzzy_low <= 0
    y_fuzzy_high = quady + y_adj - round((quady + y_adj)/y_grain) * 
        y_grain + tol * y_grain
    y_fuzzy_low = quady + y_adj - round((quady + y_adj)/y_grain) * 
        y_grain - tol * y_grain
    y_in = y_fuzzy_high >= 0 & y_fuzzy_low <= 0
    in_grid = x_in * y_in
    gridx = quadx[in_grid == 1]
    gridy = quady[in_grid == 1]
    newgridx = round((gridx + x_adj)/x_grain) * x_grain - x_adj
    newgridy = round((gridy + y_adj)/y_grain) * y_grain - y_adj
    outframe = data.frame(newgridx, newgridy, envframe[in_grid == 
        1, -c(xcol, ycol)])
    names(outframe)[1:2] = coord
    outframe
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
