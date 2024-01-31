\name{plotPath}
\alias{plotPath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotPath(fit, colors = c("gold", "green3", "blue", "brown", "pink"), logX = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{colors}{
%%     ~~Describe \code{colors} here~~
}
  \item{logX}{
%%     ~~Describe \code{logX} here~~
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
function (fit, colors = c("gold", "green3", "blue", "brown", 
    "pink"), logX = TRUE) 
{
    best.models = apply(fit$criterion.matrix, 2, which.min)
    unique.beta = unique(best.models)
    names = rep(NA, length(unique.beta))
    for (i in 1:length(unique.beta)) {
        names[i] = paste(names(best.models[best.models == unique.beta[i]]), 
            collapse = "/")
    }
    min.y = min(apply(fit$betas, 1, min))
    max.y = max(apply(fit$betas, 1, max))
    keep_models = if (logX == TRUE) 
        1:(length(fit$lambdas) - 1)
    else 1:length(fit$lambdas)
    if (logX == TRUE) {
        plot(fit$lambdas[keep_models], fit$betas[1, keep_models], 
            log = "x", type = "l", ylim = c(min.y, max.y), xlab = "LASSO penalty", 
            ylab = "Coefficients")
    }
    else {
        plot(fit$lambdas[keep_models], fit$betas[1, keep_models], 
            type = "l", ylim = c(min.y, max.y), xlab = "LASSO penalty", 
            ylab = "Coefficients")
    }
    for (i in 2:dim(fit$betas)[1]) {
        points(fit$lambdas[keep_models], fit$betas[i, keep_models], 
            type = "l")
    }
    for (i in 1:length(unique.beta)) {
        abline(v = fit$lambdas[unique.beta[i]], lwd = 3, col = colors[i])
    }
    legend("topright", names, lwd = rep(3, length(unique.beta)), 
        col = colors[1:length(unique.beta)])
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
