\name{polyNames}
\alias{polyNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
polyNames(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
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
function (X) 
{
    coef.names = dimnames(X)[[2]]
    which.poly = grep("poly\\(", coef.names)
    if (length(which.poly) > 0) {
        split1 = strsplit(coef.names[which.poly], "\\)")
        polyframe = ldply(split1, rbind)
        varframe = ldply(strsplit(unlist(strsplit(as.character(polyframe[, 
            1]), "poly\\("))[2 * (1:(length(which.poly)))], ","), 
            rbind)
        varframe = varframe[, -dim(varframe)[2]]
        varframe = data.frame(lapply(varframe, as.character), 
            stringsAsFactors = FALSE)
        if (dim(varframe)[1] == 1) {
            varframe = t(varframe)
        }
        expframe = ldply(strsplit(as.character(polyframe[, 2]), 
            "\\."), rbind)
        expframe = data.frame(lapply(expframe, as.character), 
            stringsAsFactors = FALSE)
        expframe[is.na(expframe)] = 0
        vframe = varframe[, 1:dim(expframe)[2]]
        nameframe = matrix(paste(as.matrix(vframe), "_EXP", as.matrix(expframe), 
            sep = ""), nrow(vframe), ncol(vframe))
        nameframe[expframe == "0"] = ""
        nameframe[expframe == "1"] = as.character(vframe[expframe == 
            1])
        nameframe = gsub(" ", "", nameframe)
        nameframe = as.data.frame(nameframe)
        if (is.null(dim(nameframe)) == FALSE) {
            names.out = apply(nameframe, 1, function(row) paste(row[nzchar(row)], 
                collapse = "_x_"))
        }
        if (is.null(dim(nameframe)) == TRUE) {
            names.out = nameframe
        }
        id.out = which.poly
    }
    else {
        names.out = NULL
        id.out = NULL
    }
    return(list(names = names.out, ids = id.out))
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
