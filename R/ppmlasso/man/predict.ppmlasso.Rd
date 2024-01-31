\name{predict.ppmlasso}
\alias{predict.ppmlasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.ppmlasso(object, ..., newdata, interactions = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{interactions}{
%%     ~~Describe \code{interactions} here~~
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
function (object, ..., newdata, interactions = NA) 
{
    if (any(lapply(newdata, class) == "factor")) {
        unpacknewdata = catConvert(newdata)
        newdata = unpacknewdata$X
        cat.names = setdiff(unique(unpacknewdata$cat.names), 
            NA)
        use.form = as.character(object$formula)[2]
        for (i in 1:length(cat.names)) {
            use.form = gsub(cat.names[i], paste(names(newdata)[which(unpacknewdata$cat.names == 
                cat.names[i])], collapse = " + "), use.form)
        }
        object$formula = as.formula(paste("~", use.form))
    }
    var.0 = which(apply(newdata, 2, var) == 0)
    mf = model.frame(object$formula, data = newdata)
    mt = attr(mf, "terms")
    X.des = if (!is.empty.model(mt)) 
        model.matrix(mt, mf)
    else matrix(, length(object$mu), 0L)
    X.var = X.des[, -1]
    if (is.null(object$s.means) == FALSE) {
        X.var = scale(X.var, center = object$s.means, scale = object$s.sds)
        X.des = cbind(1, X.var)
    }
    if (object$family == "area.inter") {
        if (is.na(interactions) == TRUE) {
            if (is.null(object$s.means) == FALSE) {
                X.des = cbind(X.des, min(scale(object$pt.interactions)))
            }
            if (is.null(object$s.means) == TRUE) {
                X.des = cbind(X.des, 0)
            }
        }
        if (is.na(interactions) == FALSE) {
            if (is.null(object$s.means) == FALSE) {
                X.des = cbind(X.des, scale(interactions, center = mean(object$pt.interactions), 
                  scale = sd(object$pt.interactions)))
            }
            if (is.null(object$s.means) == TRUE) {
                X.des = cbind(X.des, interactions)
            }
        }
    }
    pred.int = exp(as.matrix(X.des) \%*\% object$beta)
    return(pred.int)
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
