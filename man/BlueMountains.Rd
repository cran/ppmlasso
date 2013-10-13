\name{BlueMountains}
\docType{data}
\alias{BlueMountains}
\alias{eucalypt}
\alias{env}
\alias{availability}
\title{Blue Mountains eucalypt and environmental data.}
\description{This data set contains the observed presence locations of a Sydney eucalypt (\code{eucalypt}), the values of four
environmental variables and two variables related to site accessibility throughout the region at a spatial
resolution of 500m (\code{env}), and a matrix indicating whether locations in the region are available to the species (\code{availability}).}
\usage{data(BlueMountains)}
\format{A list with three objects:

\describe{
	\item{eucalypt}{A data frame with a column \code{X} of UTM Easting coordinates (km) and a column \code{Y} of UTM Northing coordinates (km) of observed locations of a Sydney eucalypt}
	\item{env}{A data frame containing environmental data in the Blue Mountains region near Sydney}
	\item{availability}{A 301x201 matrix with UTM Northing and Easting locations stored in \code{dimnames} indicating whether locations are accessible or not}
}

The data frame \code{env} contains the following environmental data:

\describe{
	\item{X}{UTM Easting coordinates (km)}
	\item{Y}{UTM Northing coordinates (km)}
	\item{FC}{Number of fires since 1943}
	\item{D_MAIN_RDS}{Distance from the nearest main road (m)}
	\item{D_URBAN}{Distance from the nearest urban area (m)}
	\item{RAIN_ANN}{Average annual rainfall (mm)}
	\item{TMP_MAX}{Average maximum temperature (degrees Celsius)}
	\item{TMP_MIN}{Average minimum temperature (degrees Celsius)}
}
}
\keyword{datasets}
