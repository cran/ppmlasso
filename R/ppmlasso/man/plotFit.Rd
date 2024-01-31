\name{plotFit}
\alias{plotFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotFit(fit, pred.data = data.frame(X = fit$x[fit$pres == 0], Y = fit$y[fit$pres == 0], 1, scale(fit$data[fit$pres == 0, -1], center = -fit$s.means/fit$s.sds, scale = 1/fit$s.sds)), coord = c("X", "Y"), asp = "iso", ylab = "", xlab = "", col.regions = heat.colors(1024)[900:1], cuts = length(col.regions), cex = 1.4, main.text = paste(toupper(fit$criterion), "fit"), cex.color = 1.4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{pred.data}{
%%     ~~Describe \code{pred.data} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
  \item{asp}{
%%     ~~Describe \code{asp} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{col.regions}{
%%     ~~Describe \code{col.regions} here~~
}
  \item{cuts}{
%%     ~~Describe \code{cuts} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{main.text}{
%%     ~~Describe \code{main.text} here~~
}
  \item{cex.color}{
%%     ~~Describe \code{cex.color} here~~
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
function (fit, pred.data = data.frame(X = fit$x[fit$pres == 0], 
    Y = fit$y[fit$pres == 0], 1, scale(fit$data[fit$pres == 0, 
        -1], center = -fit$s.means/fit$s.sds, scale = 1/fit$s.sds)), 
    coord = c("X", "Y"), asp = "iso", ylab = "", xlab = "", col.regions = heat.colors(1024)[900:1], 
    cuts = length(col.regions), cex = 1.4, main.text = paste(toupper(fit$criterion), 
        "fit"), cex.color = 1.4) 
{
    x.col = which(names(pred.data) == coord[1])
    y.col = which(names(pred.data) == coord[2])
    pred.int = predict(fit, newdata = pred.data[, -c(x.col, y.col)])
    levelplot(pred.int ~ pred.data[, x.col] + pred.data[, y.col], 
        asp = asp, ylab = ylab, xlab = xlab, col.regions = col.regions, 
        cuts = cuts, main = list(main.text, cex = cex), scales = list(y = list(draw = FALSE), 
            x = list(draw = FALSE)), cex = cex, colorkey = list(labels = list(cex = cex.color)))
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
