% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projective.R
\name{plus}
\alias{plus}
\alias{mult}
\alias{mult,line,numeric-method}
\alias{mult,line,matrix-method}
\alias{mult,matrix,point-method}
\title{S4 methods to add a point or multiply objects by a scalar}
\usage{
plus(x, y, ...)

mult(x, y, ...)

\S4method{mult}{line,numeric}(x, y, ...)

\S4method{mult}{line,matrix}(x, y, ...)

\S4method{mult}{matrix,point}(x, y, ...)
}
\arguments{
\item{x, y}{lines, points, gells or matrices to be added to
each other or to be displaced by a point}
}
\description{
S4 methods to add a point or multiply objects by a scalar
}
\details{
matrix operates on the left on points in line as a column vector
}
\section{Methods (by class)}{
\itemize{
\item \code{mult(x = line, y = numeric)}: multiply line by a scalar

\item \code{mult(x = line, y = matrix)}: matrix operates on line as if multiplying the points in the line on the left

\item \code{mult(x = matrix, y = point)}: matrix operates on line as if multiplying the points in the line on the left

}}
\section{Functions}{
\itemize{
\item \code{mult()}: multiply object by a scalar

}}
