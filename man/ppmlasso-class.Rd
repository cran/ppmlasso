\name{ppmlasso-class}
\Rdversion{1.1}
\docType{class}
\alias{ppmlasso-class}

\title{Class \code{"ppmlasso"}}
\description{
A class \code{ppmlasso} which represents a point process model with a LASSO penalty.
}
\section{Methods}{
  \describe{
    \item{diagnose}{\code{signature(object = "ppmlasso")}: Produce diagnostic plots for a fitted point process model.}
    \item{envelope}{\code{signature(Y = "ppmlasso")}: Produce a Monte Carlo simulation envelope for a summary function of a fitted point process model.}
    \item{predict}{\code{signature(object = "ppmlasso")}: Calculate the predicted intensity for a fitted point process model to a set of data.}
    \item{print}{\code{signature(x = "ppmlasso")}: Print the details of a fitted point process model.}
	 }
}

\author{
Ian W. Renner
}
\examples{
showClass("ppmlasso")
}
\keyword{classes}
