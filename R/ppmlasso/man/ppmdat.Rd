\name{ppmdat}
\alias{ppmdat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ppmdat(sp.xy, sp.scale, back.xy, coord = c("X", "Y"), sp.dat = getEnvVar(sp.xy = sp.xy, env.scale = sp.scale, env.grid = back.xy, coord = coord, writefile = writefile), sp.file = NA, quad.file = NA, datfilename = "PPMDat", writefile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.xy}{
%%     ~~Describe \code{sp.xy} here~~
}
  \item{sp.scale}{
%%     ~~Describe \code{sp.scale} here~~
}
  \item{back.xy}{
%%     ~~Describe \code{back.xy} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
  \item{sp.dat}{
%%     ~~Describe \code{sp.dat} here~~
}
  \item{sp.file}{
%%     ~~Describe \code{sp.file} here~~
}
  \item{quad.file}{
%%     ~~Describe \code{quad.file} here~~
}
  \item{datfilename}{
%%     ~~Describe \code{datfilename} here~~
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
function (sp.xy, sp.scale, back.xy, coord = c("X", "Y"), sp.dat = getEnvVar(sp.xy = sp.xy, 
    env.scale = sp.scale, env.grid = back.xy, coord = coord, 
    writefile = writefile), sp.file = NA, quad.file = NA, datfilename = "PPMDat", 
    writefile = TRUE) 
{
    convert = FALSE
    if (is(sp.dat, "list")) {
        convert = TRUE
        cat.names = sp.dat$cat.names
        sp.dat = sp.dat$X
    }
    if (is.character(sp.xy) == TRUE) {
        sp.file = paste(sp.xy, "Env.RData")
        load(sp.file)
    }
    if (is.na(sp.file) == FALSE) {
        load(sp.file)
    }
    if (is.na(quad.file) == TRUE) {
        dat.quad = sampleQuad(env.grid = back.xy, sp.scale = sp.scale, 
            coord = coord, writefile = writefile)
    }
    if (is.na(quad.file) != TRUE) {
        load(paste(quad.file, sp.scale, ".RData", sep = ""))
    }
    dat.quad$Pres = 0
    sp.dat$Pres = 1
    quad.x.col = which(names(dat.quad) == coord[1])
    quad.y.col = which(names(dat.quad) == coord[2])
    sp.x.col = which(names(sp.dat) == coord[1])
    sp.y.col = which(names(sp.dat) == coord[2])
    sp.var = setdiff(1:dim(sp.dat)[2], c(sp.x.col, sp.y.col))
    sp.data = data.frame(sp.dat[, sp.x.col], sp.dat[, sp.y.col], 
        sp.dat[, sp.var])
    names(sp.data) = c("X", "Y", names(sp.dat)[sp.var])
    if (any(lapply(dat.quad, class) == "factor")) {
        dat.quad = catConvert(dat.quad)$X
    }
    add.var = setdiff(names(sp.dat), names(dat.quad))
    if (length(add.var) > 0) {
        quad.names = names(dat.quad)
        for (add in 1:length(add.var)) {
            dat.quad = data.frame(dat.quad, 0)
        }
        names(dat.quad) = c(quad.names, add.var)
        dat.quad = dat.quad[, match(names(dat.quad), names(sp.dat))]
    }
    quad.var = setdiff(1:dim(dat.quad)[2], c(quad.x.col, quad.y.col))
    quad.dat = data.frame(dat.quad[, quad.x.col], dat.quad[, 
        quad.y.col], dat.quad[, quad.var])
    names(quad.dat) = c("X", "Y", names(dat.quad)[quad.var])
    quad.dat = quad.dat[, match(names(sp.dat), names(quad.dat))]
    dat.ppm = rbind(sp.data, quad.dat)
    dat.ppm$wt = weights(sp.data, quad.dat, coord)
    dimnames(dat.ppm)[[1]] = 1:dim(dat.ppm)[1]
    if (convert == TRUE) {
        dat.ppm = list(X = dat.ppm, cat.names = cat.names)
    }
    if (is.na(datfilename) == FALSE) {
        save.name = paste(datfilename, ".RData", sep = "")
        if (writefile == TRUE) {
            save(dat.ppm, file = save.name)
            print(paste("Output saved in the file", save.name))
        }
    }
    dat.ppm
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
