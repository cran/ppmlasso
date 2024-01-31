\name{print.ppmlasso}
\alias{print.ppmlasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.ppmlasso(x, ..., output = c("all", "path", "model", "interaction"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
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
function (x, ..., output = c("all", "path", "model", "interaction")) 
{
    if ("all" \%in\% output) {
        if (x$family == "poisson") {
            output = c("path", "model")
        }
        if (x$family == "area.inter") {
            output = c("path", "model", "interaction")
        }
    }
    model.grammar = if (length(x$lambdas) > 1) 
        "models"
    else "model"
    family.name = if (x$family == "poisson") 
        paste("Poisson point process", model.grammar)
    else paste("Area-interaction", model.grammar, "(r =", x$r, 
        ")")
    single.model = if (x$family == "poisson") 
        "Poisson point process model"
    else paste("Area-interaction model (r =", x$r, ")")
    pen.type = "LASSO"
    if (x$gamma > 0) {
        pen.type = paste("Adaptive LASSO (gamma =", x$gamma, 
            ")")
    }
    if (x$alpha < 1) {
        pen.type = paste("Elastic Net (alpha =", x$alpha, ")")
    }
    if ("path" \%in\% output) {
        cat(paste("Regularisation path of", length(x$lambdas), 
            family.name), "\n")
        cat(paste(pen.type, "Penalties:\n"))
        print(x$lambdas)
        cat("\nCoefficients:\n")
        print(x$betas)
        cat("\nLikelihoods:\n")
        print(x$likelihoods)
        cat("\nPenalised Likelihoods:\n")
        print(x$pen.likelihoods)
    }
    if ("model" \%in\% output) {
        cat(paste("Optimal", single.model, "chosen by", x$criterion, 
            ":\n"))
        cat(paste(pen.type, "Penalty:\n"))
        print(x$lambda)
        cat("\nCoefficients:\n")
        print(x$beta)
        cat("\nLikelihood:\n")
        print(x$likelihood)
    }
    if ("interaction" \%in\% output) {
        cat(paste("Point interactions of radius", x$r, ":\n"))
        print(x$pt.interactions)
    }
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
