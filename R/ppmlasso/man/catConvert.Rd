\name{catConvert}
\alias{catConvert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
catConvert(env.frame)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env.frame}{
%%     ~~Describe \code{env.frame} here~~
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
function (env.frame) 
{
    classes = lapply(env.frame, class)
    is.cat = which(classes == "factor")
    cat.names = names(env.frame)[is.cat]
    frame.out = env.frame[, which(classes != "factor")]
    cat.out = rep(NA, dim(frame.out)[2])
    for (v in 1:length(is.cat)) {
        current.dim = dim(frame.out)[2]
        factors = sort(unique(env.frame[, is.cat[v]]))
        fac.frame = c()
        for (i in 1:length(factors)) {
            fac.frame = cbind(fac.frame, as.numeric(env.frame[, 
                is.cat[v]] == factors[i]))
        }
        fac.frame = data.frame(fac.frame)
        frame.out = data.frame(frame.out, fac.frame)
        names(frame.out)[(current.dim + 1):dim(frame.out)[2]] = paste(cat.names[v], 
            ".", factors, sep = "")
        cat.out = c(cat.out, rep(cat.names[v], dim(fac.frame)[2]))
    }
    return(list(X = frame.out, cat.names = cat.out))
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
