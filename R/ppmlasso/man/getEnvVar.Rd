\name{getEnvVar}
\alias{getEnvVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getEnvVar(sp.xy, env.grid, env.scale, coord = c("X", "Y"), envfilename = "SpEnvData", tol = 0.01, writefile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.xy}{
%%     ~~Describe \code{sp.xy} here~~
}
  \item{env.grid}{
%%     ~~Describe \code{env.grid} here~~
}
  \item{env.scale}{
%%     ~~Describe \code{env.scale} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
  \item{envfilename}{
%%     ~~Describe \code{envfilename} here~~
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
function (sp.xy, env.grid, env.scale, coord = c("X", "Y"), envfilename = "SpEnvData", 
    tol = 0.01, writefile = TRUE) 
{
    env.grid = griddify(env.grid, tol = tol, coord = coord)
    convert = FALSE
    if (any(lapply(env.grid, class) == "factor")) {
        convert = TRUE
        out.grid = catConvert(env.grid)
        env.grid = out.grid$X
    }
    x.dat = sp.xy[, which(names(sp.xy) == coord[1])]
    y.dat = sp.xy[, which(names(sp.xy) == coord[2])]
    x.back = env.grid[, which(names(env.grid) == coord[1])]
    y.back = env.grid[, which(names(env.grid) == coord[2])]
    x.col = which(names(env.grid) == coord[1])
    y.col = which(names(env.grid) == coord[2])
    var.col = setdiff(1:dim(env.grid)[2], c(x.col, y.col))
    s.res = spatRes(env.grid)
    sp.dat = as.data.frame(matrix(NA, length(x.dat), length(var.col)))
    names(sp.dat) = names(env.grid[var.col])
    for (var in 1:length(var.col)) {
        loop.scale = min(c(s.res$x.step, s.res$y.step))
        loc = which(is.na(sp.dat[, var]))
        while (sum(is.na(sp.dat[, var])) > 0) {
            loc = which(is.na(sp.dat[, var]))
            sp.dat[loc, var] = interp(sp.xy[loc, ], loop.scale, 
                env.grid[, var.col[var]], env.grid, coord = c("X", 
                  "Y"))
            loop.scale = loop.scale * 2
        }
        cat(paste("Calculating species environmental data for variable:", 
            names(sp.dat)[var], "\n"))
        flush.console()
    }
    sp.dat = data.frame(x.dat, y.dat, sp.dat)
    names(sp.dat)[1:2] = c("X", "Y")
    if (is.na(envfilename) == FALSE) {
        save.name = paste(envfilename, ".RData", sep = "")
        if (writefile == TRUE) {
            save(sp.dat, file = save.name)
            print(paste("Output saved in the file", save.name))
        }
    }
    if (convert == TRUE) {
        sp.dat = list(X = sp.dat, cat.names = out.grid$cat.names)
    }
    sp.dat
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
