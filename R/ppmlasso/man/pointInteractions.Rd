\name{pointInteractions}
\alias{pointInteractions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pointInteractions(dat.ppm, r, availability = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat.ppm}{
%%     ~~Describe \code{dat.ppm} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{availability}{
%%     ~~Describe \code{availability} here~~
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
function (dat.ppm, r, availability = NA) 
{
    if (is(dat.ppm, "list")) {
        dat.ppm = dat.ppm$X
    }
    if (any(is.na(availability))) {
        availability = makeMask(dat.ppm)
    }
    cat(paste("Calculating point interactions", "\n"))
    flush.console()
    occupied = matrix(0, dim(availability)[1], dim(availability)[2])
    rownames(occupied) = rownames(availability)
    colnames(occupied) = colnames(availability)
    x.mat = availability
    y.mat = availability
    for (i in 1:dim(x.mat)[1]) {
        x.mat[i, ] = as.numeric(colnames(availability))
    }
    for (i in 1:dim(y.mat)[2]) {
        y.mat[, i] = as.numeric(rownames(availability))
    }
    grain = as.numeric(colnames(availability)[2]) - as.numeric(colnames(availability)[1])
    pres.x = dat.ppm$X[dat.ppm$Pres > 0]
    pres.y = dat.ppm$Y[dat.ppm$Pres > 0]
    quad.x = dat.ppm$X[dat.ppm$Pres == 0]
    quad.y = dat.ppm$Y[dat.ppm$Pres == 0]
    quad.int = rep(0, length(quad.x))
    for (i in 1:length(pres.x)) {
        sub.col = which(as.numeric(colnames(occupied)) >= pres.x[i] - 
            (r + grain) & as.numeric(colnames(occupied)) <= pres.x[i] + 
            (r + grain))
        sub.row = which(as.numeric(rownames(occupied)) >= pres.y[i] - 
            (r + grain) & as.numeric(rownames(occupied)) <= pres.y[i] + 
            (r + grain))
        sub.occ = occupied[sub.row, sub.col]
        sub.x = x.mat[sub.row, sub.col]
        sub.y = y.mat[sub.row, sub.col]
        sub.occ[(sub.x - pres.x[i])^2 + (sub.y - pres.y[i])^2 < 
            r^2] = sub.occ[(sub.x - pres.x[i])^2 + (sub.y - pres.y[i])^2 < 
            r^2] + 1
        occupied[sub.row, sub.col] = sub.occ
        quad.cells = which((quad.x - pres.x[i])^2 + (quad.y - 
            pres.y[i])^2 <= (2 * r)^2)
        quad.int[quad.cells] = quad.int[quad.cells] + 1
    }
    int.q = rep(0, length(quad.int))
    for (quad.i in which(quad.int > 0)) {
        sub.col = which(as.numeric(colnames(occupied)) >= quad.x[quad.i] - 
            (r + grain) & as.numeric(colnames(occupied)) <= quad.x[quad.i] + 
            (r + grain))
        sub.row = which(as.numeric(rownames(occupied)) >= quad.y[quad.i] - 
            (r + grain) & as.numeric(rownames(occupied)) <= quad.y[quad.i] + 
            (r + grain))
        sub.occ = occupied[sub.row, sub.col]
        sub.availability = availability[sub.row, sub.col]
        sub.x = x.mat[sub.row, sub.col]
        sub.y = y.mat[sub.row, sub.col]
        sub.cell = (sub.x - quad.x[quad.i])^2 + (sub.y - quad.y[quad.i])^2 <= 
            r^2 & sub.availability > 0
        int.q[quad.i] = sum(sub.occ[sub.cell] > 0, na.rm = TRUE)/sum(sub.cell, 
            na.rm = TRUE)
    }
    int.p = rep(0, length(pres.x))
    for (pres.i in 1:length(pres.x)) {
        sub.col = which(as.numeric(colnames(occupied)) >= pres.x[pres.i] - 
            (r + grain) & as.numeric(colnames(occupied)) <= pres.x[pres.i] + 
            (r + grain))
        sub.row = which(as.numeric(rownames(occupied)) >= pres.y[pres.i] - 
            (r + grain) & as.numeric(rownames(occupied)) <= pres.y[pres.i] + 
            (r + grain))
        sub.occ = occupied[sub.row, sub.col]
        sub.availability = availability[sub.row, sub.col]
        sub.x = x.mat[sub.row, sub.col]
        sub.y = y.mat[sub.row, sub.col]
        sub.cell = (sub.x - pres.x[pres.i])^2 + (sub.y - pres.y[pres.i])^2 <= 
            r^2 & sub.availability > 0
        int.p[pres.i] = sum(sub.occ[sub.cell] > 1, na.rm = TRUE)/sum(sub.cell, 
            na.rm = TRUE)
    }
    interactions = c(int.p, int.q)
    interactions
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
