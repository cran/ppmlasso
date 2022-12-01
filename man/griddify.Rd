\name{griddify}
\alias{griddify}
\title{Ensure that a geo-referenced matrix of environmental grids is rectangular}
\description{This function ensures that the coordinates of the supplied geo-referenced matrix of environmental grids constitutes a rectangular grid.}
\usage{
griddify(envframe, tol = 0.01, coord = c("X", "Y"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{envframe}{The geo-referenced matrix of environmental grids.}
  \item{tol}{The tolerance level within which to correct coordinate errors.}
  \item{coord}{A vector containing the names of the longitude and latitude coordinates.}
}
\details{The functions in the \code{ppmlasso} package require a set of quadrature points along a rectangular grid. At times a set of quadrature points with a desired spatial resolution of \eqn{x_\delta \times y_\delta} will have some minor machine error in some coordinates such that the coordinates as supplied do not consistute a rectangular grid. The \code{griddify} function corrects this error as follows:

Let \eqn{\{x_1, x_2, \ldots, x_n\}} and \eqn{\{y_1, y_2, \ldots, y_n\}} be the supplied coordinates contained in \code{envframe}. The function first determines the spatial resolution \eqn{x_\delta \times y_\delta} based on the median of the differences in the unique values of \eqn{x_i} and \eqn{y_i} as well as the coordinates of a rectangular grid with this spatial resolution \eqn{\{x^{grid}_1, x^{grid}_2, \ldots, x^{grid}_n\}} and \eqn{\{y^{grid}_1, y^{grid}_2, \ldots, y^{grid}_n\}}. Given the tolerance \eqn{\epsilon} supplied to \code{tol}, any coordinate \eqn{x_i} for which \eqn{0 < \left|x_i - x^{grid}_i\right| \leq \epsilon \times x_\delta} will be adjusted to \eqn{x^{grid}_i}. Likewise, any coordinate \eqn{y_i} for which \eqn{0 < \left|y_i - y^{grid}_i\right| \leq \epsilon \times y_\delta} will be adjusted to \eqn{y^{grid}_i}.

Any environmental variables contained in \code{envframe} are left unchanged.
}
\value{A data frame containing the coordinates on a rectangular grid as well as any environmental variables left unchanged.}
\author{Ian W. Renner}
\examples{
X = seq(0, 5, 1)
Y = seq(1, 11, 2)
XY = expand.grid(X, Y) # generate 1 x 2 rectangular grid
names(XY) = c("X", "Y")
#move some coordinates off of rectangular grid
XY$X[1] = XY$X[1] - 0.01
XY$Y[1] = XY$Y[1] - 0.01
XY$X[7] = XY$X[7] + 0.01
XY$Y[7] = XY$Y[7] + 0.01

#generate environmental variables
XY$V1 = 0.1*XY$X + 0.2*XY$Y + rnorm(36, 0, 1)
XY$V2 = -0.2*XY$X + 0.1*XY$Y + 0.05*XY$X*XY$Y + rnorm(36, 0, 5)

XY_grid = griddify(XY)
}