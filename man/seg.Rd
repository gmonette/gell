% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projective.R
\name{seg}
\alias{seg}
\title{Create Segments}
\usage{
seg(x, ...)
}
\arguments{
\item{x, ...}{lines, points or \code{gell}s to be combined in a matrix of points
with class \code{segment}}
}
\description{
Segments are constructed from sequences of points that are intended to
be drawn a lines joining the points.
}
\details{
A single segment, defined by two points, is terminated at the the point,
in contrast with \code{line}s that are not terminated. Segments can be terminated
at a point at infinity and, in that case, will not appear terminated when
plotted. A sequence of segments can include a point or points are infinity.
}
\examples{
     
  l(rbind(c(1,1,1), l(1,-1,1)))
  join(l(1,1,1), l(1,-1,1)) \%>\% 
    rbind(p(2,3,1)) \%>\% p \%>\%  seg
}
