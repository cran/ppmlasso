\name{sampleQuad}
\alias{sampleQuad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleQuad(env.grid, sp.scale, coord = c("X", "Y"), quadfilename = "Quad", tol = 0.01, writefile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env.grid}{
%%     ~~Describe \code{env.grid} here~~
}
  \item{sp.scale}{
%%     ~~Describe \code{sp.scale} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
  \item{quadfilename}{
%%     ~~Describe \code{quadfilename} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{writefile}{
%%     ~~Describe \code{writefile} here~~
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
function (env.grid, sp.scale, coord = c("X", "Y"), quadfilename = "Quad", 
    tol = 0.01, writefile = TRUE) 
{
    convert = FALSE
    env.grid = griddify(env.grid, tol = tol, coord = coord)
    if (any(lapply(env.grid, class) == "factor")) {
        convert = TRUE
        out.grid = catConvert(env.grid)
        env.grid = out.grid$X
    }
    x.col = which(names(env.grid) == coord[1])
    y.col = which(names(env.grid) == coord[2])
    res = spatRes(env.grid, coord = coord)
    x.step = res$x.step
    y.step = res$y.step
    f.name = list()
    quad.list = list()
    for (i in 1:length(sp.scale)) {
        i.scale = sp.scale[i]
        x.o = min(env.grid[, x.col]) - floor(min(env.grid[, x.col])/x.step) * 
            x.step
        y.o = min(env.grid[, y.col]) - floor(min(env.grid[, y.col])/y.step) * 
            y.step
        if (x.o/x.step > 0.5) {
            x.o = x.o - x.step
        }
        if (y.o/y.step > 0.5) {
            y.o = y.o - y.step
        }
        is.on.scale = abs((env.grid[, x.col]/i.scale) - round(env.grid[, 
            x.col]/i.scale) - x.o/i.scale) + abs((env.grid[, 
            y.col]/i.scale) - round(env.grid[, y.col]/i.scale) - 
            y.o/i.scale) < 1e-08
        dat.quad = env.grid[is.on.scale, ]
        dec.x = max(unlist(lapply(dat.quad[, x.col], decimalCount)))
        dec.y = max(unlist(lapply(dat.quad[, y.col], decimalCount)))
        dat.quad[, x.col] = unlist(lapply(dat.quad[, x.col], 
            zapCoord, dec.x))
        dat.quad[, y.col] = unlist(lapply(dat.quad[, y.col], 
            zapCoord, dec.y))
        if (is.na(quadfilename) == FALSE) {
            f.name[[i]] = paste(quadfilename, sp.scale[i], ".RData", 
                sep = "")
            if (writefile == TRUE) {
                save(dat.quad, file = f.name[[i]])
                print(paste("Output saved in the file", f.name[[i]]))
            }
        }
        if (length(sp.scale) > 1) {
            quad.list[[i]] = dat.quad
        }
    }
    if (convert == TRUE) {
        dat.quad = list(X = dat.quad, cat.names = out.grid$cat.names)
    }
    if (length(sp.scale) == 1) 
        dat.quad
    else quad.list
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
