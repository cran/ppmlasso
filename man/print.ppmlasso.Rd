\name{print.ppmlasso}
\alias{print.ppmlasso}
\title{Print a fitted regularisation path}
\description{This function prints output from a fitted regularisation path.}
\usage{
\method{print}{ppmlasso}(x, \dots, output = c("all", "path", "model", "interaction"))
}
\arguments{
  \item{x}{A regularisation path fitted by \code{\link{ppmlasso}}.}
  \item{\dots}{Further arguments controlling the printed output.}
  \item{output}{This argument controls what output is printed to the screen. If \code{output} includes \code{"path"},
information about the entire regularisation path is printed. If \code{output} includes \code{"model"},
information about the model that optimises the given criterion is printed. If \code{output} includes \code{"interaction"},
information about the point interactions is printed. Setting \code{output = "all"} will print all available information.}
}
\value{N/A}
\author{Ian W. Renner}
\seealso{
\code{\link{ppmlasso}} for fitting regularisation paths.
}
\examples{
# Fit a regularisation path of Poisson point process models
data(BlueMountains)
ppm.form = ~ poly(FC, TMP_MIN, TMP_MAX, RAIN_ANN, degree = 2)
ppm.fit  = ppmlasso(ppm.form, sp.xy = BlueMountains$eucalypt, 
env.grid = BlueMountains$env, sp.scale = 1, n.fits = 20, writefile = FALSE)
print(ppm.fit)
}
