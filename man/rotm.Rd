% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projective.R
\name{rotm}
\alias{rotm}
\title{Rotation matrix}
\usage{
rotm(theta = pi/2)
}
\arguments{
\item{theta}{angle in radians that matrix would rotate}
}
\description{
Rotation matrix
}
\examples{
 library(spida2)
 rotm(pi/100)
 init()
 p(0,0,1) \%>\% pl
 l(1,1,1) \%>\% pl
 l(1,1,1) \%>\% mult(rotm(.1)) \%>\% pl
 l(1,1,1) \%>\% mult(rotm(pi/3)) \%>\% pl
 p(2,1,1) \%>\% {mult(rotm(pi/3),.)} \%>\% pl
 p(2,1,1) \%>\% {mult(rotm(0/3),.)} \%>\% pl
 
 
 l(1,1,0) \%>\% pl
 l(1,1,0) \%>\% mult(cbind(c(1,0),c(2,1))) \%>\% pl
 l(1,1,0) \%>\% mult(cbind(c(1,0),c(3,1))) \%>\% pl
 p(0,0,1) \%>\% pl(pch = 16)
 l(1,1,1) \%>\% pl
 l(1,1,1) \%>\% mult(cbind(c(1,0),c(3,1))) \%>\% pl
 l(1,1,1) \%>\% mult(cbind(c(0,-1),c(1,0))) \%>\% pl
}
