\name{findRes}
\alias{findRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findRes(scales, lambda = 0, coord = c("X", "Y"), sp.xy, env.grid, formula, tol = 0.01, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scales}{
%%     ~~Describe \code{scales} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
  \item{sp.xy}{
%%     ~~Describe \code{sp.xy} here~~
}
  \item{env.grid}{
%%     ~~Describe \code{env.grid} here~~
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (scales, lambda = 0, coord = c("X", "Y"), sp.xy, env.grid, 
    formula, tol = 0.01, ...) 
{
    form.1 = formula
    likelihoods = rep(NA, length(scales))
    sp.data = getEnvVar(sp.xy = sp.xy, env.scale = min(scales), 
        env.grid = env.grid, coord = coord, envfilename = NA, 
        tol = tol)
    for (sc in 1:length(scales)) {
        formula = form.1
        if (lambda == 0) {
            data = ppmdat(sp.xy = sp.xy, sp.scale = scales[sc], 
                back.xy = env.grid, sp.dat = sp.data, sp.file = NA, 
                quad.file = NA, writefile = FALSE)
            if (is(data, "list")) {
                use.form = as.character(formula)[2]
                cat.names = setdiff(unique(data$cat.names), NA)
                for (i in 1:length(cat.names)) {
                  use.form = gsub(cat.names[i], paste(names(data$X)[which(data$cat.names == 
                    cat.names[i])], collapse = " + "), use.form)
                }
                formula = as.formula(paste("~", use.form))
                data = data$X
            }
            glm.form = as.formula(paste("Pres/wt ~ ", as.character(formula)[2], 
                sep = ""))
            glm.fit = suppressWarnings(glm(glm.form, data = data, 
                weights = data$wt, family = poisson()))
            eps = 1e-09
            while (glm.fit$deviance > glm.fit$null.dev) {
                glm.fit = suppressWarnings(glm(glm.form, data = data, 
                  weights = data$wt, family = poisson(), control = list(epsilon = eps)))
                eps = eps/10
            }
            likelihoods[sc] = sum(data$wt * (data$Pres/data$wt * 
                log(glm.fit$fitted) - glm.fit$fitted)) - sum(log(1:sum(data$Pres > 
                0)))
        }
        if (lambda != 0) {
            sc.fit = ppmlasso(glm.form, sp.scale = scales[sc], 
                lamb = lambda, data = ppmdat(sp.xy = sp.xy, sp.scale = scales[sc], 
                  back.xy = env.grid, sp.dat = sp.data, sp.file = NA, 
                  quad.file = NA, writefile = FALSE), ...)
            likelihoods[sc] = sc.fit$pen.likelihood[1]
        }
    }
    plot(scales, likelihoods, log = "x", type = "o", pch = 16, 
        xlab = "Spatial Resolution", ylab = "Likelihood")
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
