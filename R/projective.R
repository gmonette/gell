#' gell: Working with Ellipses
#' 
#' This is a part of a package for plotting objects
#' represented in projective 2 dimensional space
#' Using objects 'point' and 'line'
#' 
#' General approach
#' 
#' See summary and examples at end
#' 
#' Plotting Functions:
#' 
#' - [pl()]  a method that plots most objects
#' - [text()] a method that adds text at a location
#' - [gell::text] at a point
#' - [text(point)] at a point
#' - [gell::contour.function] S3 method using [graphics::contour]
#' 
#' Methods:
#' 
#' - fta: 'from', 'to', 'along'
#' - axes: gell, point (dir), type
#' 
#' Idea:
#' 
#' - Use the projective representation of points and lines
#'   to include the points at infinity and the line at infinity
#'   thereby making it easy to work with parallel lines and, e.g.
#'   the line joining a point with the intersection of two lines
#'   that happen to be parallel, etc.
#' - Manipulated projective points and lines with a suite of
#'   operators that coerce finite points and lines to projective
#'   form.
#' - We need three classes of objects (but class might be implicit)
#'   - ppoints
#'   - plines (unbounded)
#'   - psegs  for segments and, e.g. parallelograms, etc.
#' - A sequence of n points can be intererpreted as n-1 segments
#'   or as n-1 lines
#' - A sequence of n lines can generate n-1 intersections
#' - Generalized ellipses, class 'gell'
#'   For now represented using a svd representation using
#'   orthogonal \eqn{$\Gamma$} and singular values which may be
#'   0 or Inf. 
#' - Write plotting methods for base R and for lattice that
#'   plot projective points and lines and gell's. Note that lines can be
#'   drawn with 'abline' to avoid problems with endpoints.
#'  
#'    
#' 
#' - Use projective geometry to easily compute the intersection of two lines
#'   and the line joining two points
#' - Issues:
#'   - Lines can be represented as a projective line (3-vector) or as
#'     two endpoints.  The former needs to be transformed to the latter
#'     for actual plotting. 
#'   - We need to distinguish between a line segment and a line.
#'     - lines could be plotted with `abline` so no endpoints
#'     - line segments can be turned into lines
#'  
#' Operators:
#' 
#' - 1 or more lines
#'   - intersection
#'   - copuntal
#'   - rotate
#'   - transform
#'   - 
#'   - dual point of a line (o a finite line = point at inf perpendicular)
#'
#' - line and point
#'   - perpendicular line through point (might add ell to get conjugate line)
#'   - Relation: point in line?
#' - 2 or more points
#'   - sum
#'   - scalar
#' - [gell()] and
#'   - point
#'     - scale to gell
#'     - conjugate [axes()]
#'     - conjugate parallelogram
#'     - subtending rectangle
#' - gell and
#'   - line
#'     - tangent point
#' - gell and gell
#'   - locus of osculation
#'        
#' @docType package    
#' @name gell
#' @keywords internal
"_PACKAGE"

## FUNCTIONS ####
library(spida2)
library(magrittr)
setOldClass('gell')
setOldClass('point')
setOldClass('line')
setOldClass('ell')
setOldClass('seg')
# library(pracma)
# library(spida2)
# makeActiveBinding("ans", function() .Last.value, .GlobalEnv)
# who <- pracma::who
# ?makeActiveBinding
# COLS <- c('red','cyan','purple','yellow','teal' ,'blue','#4444FF','green','goldenrod',
#                '#23423411',
#                '#014488')
# red <- 'red'
#   blue <- 'blue'
#     green <- 'green'
#       black <- 'black'
# 
#' Initialize plotting surface
#' 
#' @param from lower left corner
#' @param to upper right corned
#' @param fac expand c(-1,1)
#' @param xlab x label
#' @param ylab y label
#' @param asp (default: 1)
#' @export       
init <- function(xlim = fac*c(-1,1), ylim = xlim, fac = 3,
                 xlab = '', ylab = '', asp = 1,...) {
  plot(0,0, xlim = xlim, ylim = ylim, type ='n', xlab = xlab, ylab = ylab, asp = asp,...)
  abline(h=0, col = "#44444422")
  abline(v=0, col = "#44444422")
}
#' Append to a List
#' 
#' Needs more work
#' 
#' @param l a list
#' @param ... objects to be appended
#' @export
clist <- function(l, ...) {
  # add to spida2
  a <- list(...)
  if(is.null(names(a))) names(a) <-  rep('', length(a))
  if(is.null(names(l))) nl <- rep('', length(l))
  nl <- names(l)
  first <- length(l) + 1
  last <- length(l) + length(a)
  l[first:last] <- a
    names(l) <- c(nl,names(a))
  l
}

if(F) {
  za <- list(a= 'a',b= 'b')
  zb <- list(d= 'd', e = 'e')
  zr <- list(a = 'z', f = 'r')
  
  c(za, zb)  
  c(list(za), x=list(zb), zr)  
  clist(za, z=zb)
  clist(za, z=zb, za, zr) %>% names
}
#' Flip Rows of a Matrix or Other Object
#' 
#' @param x object to be flipped
#' @returns reverses most objects, e.g. lists, reverses the rows of a matrix
#' @export
flip <- function(x,...) {
  UseMethod('flip')
}           
#' @export
flip.default <- function(x, ind = rev(seq_along(x)), ...){
  ret <- x[ind]
  class(ret) <- class(x)
  ret
}
#' @export
flip.matrix <- function(x, ind = rev(seq_len(nrow(x))),...){
  ret <- x[ind,,drop = FALSE]
  class(ret) <- class(x)
  ret
}
#' @export
flip.seg <- function(x,ind = rev(seq_len(nrow(x))),...){
  ret <- x[ind,,drop = FALSE]
  class(ret) <- class(x)
  ret
}
#' @export
flip.line <- function(x, ind = rev(seq_len(nrow(x))), ...){
  ret <- x[ind,,drop = FALSE]
  class(ret) <- class(x)
  ret
}
#' @export
flip.point <- function(x, ind = rev(seq_len(nrow(x))), ...){
  ret <- x[ind,,drop = FALSE]
  class(ret) <- class(x)
  ret
}




# p(0,0) %>% rbind(p(1,1)) %>% flip

#' Normalize homogeneous coordinates or gell
#' 
#' @param x an object represented as a n x 3 numerical matrix
#' 
#' @export
n <- function(x,...){
  UseMethod('n')
}
#' @export
n.line <- function(x) {
  if(ncol(x) != 3) stop('can only normalize homogeneous 3-vectors')
  div <- ifelse(x[,3] == 0, pracma::hypot(x[,1],x[,2]), x[,3])
  ret <- x/div
  class(ret) <- class(x)
  ret
}

#' @export
n.point <- function(x) {
  if(ncol(x) != 3) stop('can only normalize homogeneous 3-vectors')
  div <- ifelse(x[,3] == 0, pracma::hypot(x[,1],x[,2]), x[,3])
  ret <- x/div
  class(ret) <- class(x)
  ret
}

#' @export
n.matrix <- function(x){
  if(ncol(x) != 3) stop('can only normalize homogeneous 3-vectors')
  div <- ifelse(x[,3] == 0, pracma::hypot(x[,1],x[,2]), x[,3])
  ret <- x/div
  class(ret) <- class(x)
  ret
}

#' @export
n.default <- function(x) {
  if(ncol(x) != 3) stop('can only normalize homogeneous 3-vectors')
  div <- ifelse(x[,3] == 0, pracma::hypot(x[,1],x[,2]), x[,3])
  ret <- x/div
  class(ret) <- class(x)
  ret
}

#' @export
n.gell <- function(x){
  if(max(x$d) < Inf) {
    A <- x$gamma %*% diag(x$d)
    ret <- gell(A=A, radius = x$radius, center = z$center)
  } else if (min(x$d) < Inf) {
    fin <- x$d < Inf
    a1 <- x$gamma[,!fin]
    g1 <- a1/pracma::hypot(a1[1],a1[2])
    g2 <- c(g1[2], - g1[1])
    a2 <- x$gamma[,fin] * x$d[fin]
    dfin <- sum(g2 * a2)
    if(dfin < 0){
      dfin <- -dfin
      g2 <- -g2
    }
    ret <- gell(gamma = cbind(g1,g2), d = c(Inf, dfin), radius = x$radius)
  } else {
    ret <- x
  }
  ret
}

if(FALSE) {
  init()
  gs <- gell(A=cbind(c(1,1),c(0,0)))
  gs %>% pl
  gs %>% dual %>% n %>% pl
}    


#' Plot generic function

#' @export
pl <- function(x, ...) {
  # plot method
  UseMethod('pl')
}

#' @export
r <- function(...) {  # maybe get rid of this
  rbind(c(...))
}

##+ point ####
setOldClass('point')
#' Constructor for point object
#' 
#' Points are represented by rows of n x 3 matrices where
#' each row is the projective representation of a point
#' in homogeneous coordinates.
#' 
#' @details
#' The row vector (a, b, c) represents the point consisting of
#' the affine point (a/c, b/c) if c is not equal to 0, or
#' the point at infinity (a, b, 0) representing the 'direction'
#' pointing from the origin through the point (a, b).
#' 
#' Two parallel lines meet at a point at infinity
#' 
#' The vector (0, 0, 0) does not define a point. 
#' 
#' @param x,... one of the following:
#' * a vector of length 2 representing an affine point (i.e. not a point at infinity)
#' * a vector of length 3: if the third element is not zero: an affine point in homegeneous coordinates or,
#'   if the third element is 0: a point at infinity, i.e. a direction.
#' * a n by 2 matrix of affine points.
#' * a n by 3 matrix of both affine points and points at infinity
#' * objects that can be interpreted as projective points, e.g. objects of class 'seg' 
#' @returns
#' An object of class `point` with the structure of an n by 3 matrix where each row
#' represents a point.
#' 
#' @export
p <- function(...) {
  x <- list(...)
  if(length(x) %in% 2:3 && length(unlist(x)) == length(x)) {
    return(p(unlist(x)))
  }
  # disp(x)
  nn <- names(x)
  # disp(nn)
  x <- lapply(x, function(x) {
    if(is.matrix(x)) {
      if(ncol(x) == 2) {
        cbind(x, 1)
      } else {
        x
      } 
    } else if(length(x) == 2) {
      rbind(c(x,1))
    } else {
      rbind(c(x))
    }
  }
  )
  x <- lapply(seq_along(x), function(ii) {
    if(is.null(row.names(x[[ii]]))) rownames(x[[ii]])<-rep(nn[ii],nrow(x[[ii]]))
    x[[ii]]
  })
  x <- do.call(rbind, x)
  x <- n(x)
  class(x) <- 'point'
  x
}

# p <- function(x, ...) {
#   UseMethod('p')
# } 
# p.default <- function(x, ...) {
#   # possibilities: x is a
#   # vector of length 2 or 3 followed by other objects of length 2 or 3
#   
#   
#   if(length(x) == 1) x <- rbind(c(x, ...))
#   else x <- rbind(x[,],...)
#   if(ncol(x) == 2) x <- cbind(x,1)
#   x <- n(x)
#   class(x) <- 'point'
#   x
# }
##+ point ####
setOldClass('point')
#' Constructor for line object
#' 
#' Lines are represented by rows of n x 3 matrices where
#' each row is the projective representation of a 2D line
#' in homogeneous coordinates, or the representation of the
#' line at infinity, (0, 0, c) where c is not zero.
#' 
#' @details
#' The row vector (a, b, c) represents the line consisting
#' of pairs (x, y) such that ax + by + c = 0.
#' 
#' Note that if c = 0, the line goes through the origin.
#' 
#' The line at infinity is represented by the vector (0, 0, c),
#' with c != 0. 
#' 
#' The vector (0, 0, 0) does not define a line. 
#' 
#' @param x,... one of the following:
#' * a vector of length 2 representing a line through the origin. The 45 degree
#'   line is represented by (1, -1, 0), i.e. the points (x, y) satisfying 1*x + (-1)*y + 0 = 0.
#' * a vector of length 3: if the third element is not zero: a line not through
#'   the origin point in homegeneous coordinates or,
#'   if the third element is 0: a point at infinity, i.e. a direction.
#' * a n by 2 matrix of affine points.
#' * a n by 3 matrix of both affine points and points at infinity
#' * objects that can be interpreted as projective points, e.g. objects of class 'seg' 
#' @returns
#' An object of class `line` with the structure of an n by 3 matrix where each row
#' represents a line.
#' 
#' @export
l <- function(x,...) {
  # constructor for 'line'
  UseMethod('l')
} 

#' @export
l.line <- function(x,...) {
  n(x)
}

#' @export
l.default <- function(x, ...) {
  if(length(x) == 1) x <- rbind(c(x, ...))  # why???
  else x <- rbind(x,...)
  if(ncol(x) == 2) x <- cbind(x,0)  # line through origin
  x <- n(x)           # normalize
  class(x) <- 'line'
  x
}

#' @export
l.seg <- function(sg) {
  if(nrow(sg) <2) stop('need at least two points')
  ret <- matrix(0, nrow = 0, ncol = 3)
  for(i in 2:nrow(sg)){
    ret <- rbind(ret, join(p(sg[i-1,]),p(sg[i,])))
  }
  l(ret)
}

#############  SEGMENTS  #################################
#' Create Segments
#' 
#' Segments are constructed from sequences of points that are intended to
#' be drawn a lines joining the points.  
#' 
#' A single segment, defined by two points, is terminated at the the point,
#' in contrast with `line`s that are not terminated. Segments can be terminated
#' at a point at infinity and, in that case, will not appear terminated when 
#' plotted. A sequence of segments can include a point or points are infinity.
#' 
#' @param x,... lines, points or `gell`s to be combined in a matrix of points
#'        with class `segment`
#' @examples      
#'   l(rbind(c(1,1,1), l(1,-1,1)))
#'   join(l(1,1,1), l(1,-1,1)) %>% 
#'     rbind(p(2,3,1)) %>% p %>%  seg
#' @export
seg <- function(x, ...) {
  UseMethod('seg')
}
#' @export
seg.line <- function(ll, ...) {
  lfrom <- l(ll)
  if(nrow(lfrom) < 2) stop('need at least two lines')
  ind <- seq_len(nrow(lfrom))
  lto <- l(rbind(ll[-1,],ll[1,]))
  pts <- join(lfrom, lto)
  class(pts) <- 'seg'
  pts
}

#' @export
seg.point <- function(pp, close = FALSE, ...) {
  pp <- p(pp)
  if(close) pp <- close.seg(pp)
  class(pp) <- 'seg'
  pp
}

#' @export
seg.default <- seg.point

#' @export
seg.gell <- function(ge, ...) {
  if(max(ge$d) < Inf) {
    U <- circle()
    pts <- t(ge$gamma %*% (t(U) * ge$radius * ge$d ) + ge$center)
    seg(pts)
  } else {
    imax <- which(ge$d == Inf)
    imin <- which(ge$d != Inf)
    pts <- rbind(
      c(ge$gamma[,imax],0),
      c(ge$gamma[,imin]*ge$radius*ge$d[imin],1),
      c(-ge$gamma[,imax],0),
      c(-ge$gamma[,imin]*ge$radius*ge$d[imin],1),
      c(ge$gamma[,imax],0)
    )
  }
  seg(pts)
}
#' Close a sequence of segments
#' 
#' Terminate a sequence of segments with the first point in the sequence
#' 
#' @param sg segments to be terminated with first point
#' @describeIn seg close sequence of points by ending with first point
#' @export
close.seg <- function(sg) {
  if(!identical(sg[1,],sg[nrow(sg),])) sg <- rbind(sg, sg[1,])
  seg(p(sg))
}

#' 
setOldClass('line')
setOldClass('point')
#' S4 Generic Function for the Join of Projective Lines and Points
#' 
#' The intersection of two lines or the line generated by two points
#' 
#' @param x,y two projective lines or two projective points
#' @returns the line joining two points or the intersection of
#'          two lines. The the lines are parallel the intersection
#'          is a point at infinity.
#' @examples
#' join(p(c(1,1),c(1,-1)), p(c(2,2)))
#' join(
#'   p(1,1,0),  # point at infinity)
#'   p(c(1,2))
#' )
#' join(p(c(1,1,1)), p(c(1,1,1))) %>% join(l(c(1,-1,0))) %>% join(p(c(1,1,1)))
#' @export
setGeneric('join', function(x, y, ...) standardGeneric('join'))

## join: point point -------------------
#' @export
setMethod('join', signature('point','point'),
          function(x, y, ...){
            # recycle o2 if fewer rows
            a1 <- x
            a2 <- y
            n1 <-nrow(a1)
            n2 <-nrow(a2)
            ind2 <- rep(seq_len(n2),length.out=n1)
            ret <- pracma::cross(a1,a2[ind2,, drop = FALSE])
            ret <- n(ret)
            # if('point' %in% class(o1)) class(ret) <- 'line'
            # if('line' %in% class(o1)) class(ret) <- 'point'
            # if(!missing(class.out)) class(ret) <- class.out
            l(ret)
          })

## join: line line ------------------------
#' @export
setMethod('join', signature('line','line'),
          function(x, y, ...){
            # recycle o2 if fewer rows
            a1 <- x
            a2 <- y
            n1 <-nrow(a1)
            n2 <-nrow(a2)
            
            ind2 <- rep(seq_len(n2),length.out=n1)
            
            # disp(a1)
            # disp(a2[ind2,, drop = FALSE])
            ret <- pracma::cross(a1,a2[ind2,, drop = FALSE])
            ret <- n(ret)
            # if('point' %in% class(o1)) class(ret) <- 'line'
            # if('line' %in% class(o1)) class(ret) <- 'point'
            # if(!missing(class.out)) class(ret) <- class.out
            p(ret)
          })

# test

#' S3 Plotting Generic Function for Graphic Objects in the Gell Package
#' 
#' @param x object to be plotted
#' @param ... plotting parameters
#' 
#' @export
pl <- function(x, ...) {
  UseMethod('pl')
}
#' Plot segments
#' 
#' @param close if TRUE repeat first point at end of sequence (default: FALSE)
#' @param type 'p', 'b' or 'l' (default)
#' @param fill fill color (default: NULL)
#' @param ... other parameters passed to plotting function
#' @export
pl.seg <- function(x, close = FALSE, type = 'l', fill = NULL, verbose = 0, ...) {
  # 
  # Find a radius to include all points
  # 
  x <- p(x)
  # if(x[1,3] == 0 || x[nrow(x),3] == 0) stop('First and last points should be finite.')
  screen <- rbind(
    p(matrix(par('usr'), nrow = 2)),
    x
  )
  # need to include all points wherever centered
  rad <- 3 * max(pracma::hypot(screen[, 1], screen[,2]))
  if(verbose > 0) disp(rad)
  
  # ge <- gell(radius = rad)  # for points at infinity
  
  if(close) x <- close.seg(x)
  
  if(x[1,3] == 0) {
    ret <- matrix(0, ncol = 3, nrow =0)
  } else {
    ret <- x[1,,drop = FALSE]
  }
  att <- list()     
  for(i in 2:nrow(x)) {
    if(verbose > 0) disp(i)
    # possibilities
    # - current point finite: 
    #   - previous point finite: add point to ret
    #   - previous point infinite: 
    #     find point on circle and add point on circle and current point to ret
    # - current point infinite:
    #   - previous point finite: add point on circle
    #   - previous point infinite: skip
    # TODO LATER: join points on circle to avoid visible artefacts
    if(x[i,3] != 0) {  # current point finite
      if(x[i-1,3] != 0) {                  # previous point finite
        if(verbose > 0) disp('fin-fin')
        ret <- rbind(ret, x[i,])
        if(verbose > 0) disp(ret)
      } else {                             # previous point infinite
        if(verbose > 0) disp('infin-fin')
        ge <- gell(center=x[i,1:2], radius = rad)
        if(verbose > 0) disp(ge)
        att <- append(att,list(ge))
        infp <- p(axes(ge, p(x[i-1,]))[[1]][1,])
        att <- append(att,list(infp))
        ret <- rbind(ret, infp, ppi = x[i,])
        if(verbose > 0) disp(ret)
      }
    } else {  # current point infinite
      
      if(x[i-1,3] == 0) {  # previous point infinite
        if(verbose > 0) disp('infin-infin')
      } else {  # previous point finite
        if(verbose > 0) disp('fin-infin')
        
        ge <- gell(center=x[i-1,1:2], radius = rad)
        if(verbose > 0) disp(ge)
        att <- append(att, list(ge))
        infp <- p(axes(ge, p(x[i,]))[[1]][1,])
        att <- append(att, list(infp))
        ret <- rbind(ret, infp)
        if(verbose > 0) disp(ret)
        
      }
    }
  }
  if(is.null(fill)) {
    lines(ret[,1:2, drop = FALSE], type = type, ...)
  } else {
    polygon(ret[,1],ret[,2], col = fill, ...)
  }
  if(verbose > 0) disp(att)  
  ret <- seg(p(ret))
  attr(ret,'att') <- att
  invisible(ret)
}
#' @export
pl.list <- function(x, ...) {
  lapply(x, pl, ...)
  invisible(x)
}

#' @export
pl.default <- function(x, ...) {
  cat(paste("\\npl: no method for class:", class(x), "\\n"))
  invisible(x)
}

#' @export
pl.matrix <- function(x, type = 'p', ...) {
  points(x[,1:2, drop = FALSE], type = type, ...)
  invisible(x)
}

#' @export
pl.line <- function(x,...) {
  x <- l(x)
  for(i in seq_len(nrow(x))) {
    ll <- x[i,]
    if(ll[2] != 0) abline(a = -ll[3]/ll[2], b = -ll[1]/ll[2], ...)
    else abline(v = -ll[3]/ll[1], ...)
  }
  invisible(x)
}

#' @export
pl.point <- function(x,...) {
  x <- p(x)
  x <- x/x[,3]
  points(x[,1:2,drop=FALSE], ...)
  invisible(x)
}

## plus, mult ####
#' S4 methods to add a point or multiply objects by a scalar
#' 
#' @param x,y lines, points, gells or matrices to be added to
#'        each other or to be displaced by a point
#'        
#' @export        
setGeneric('plus', function(x,y,...) standardGeneric('plus') )
#' @describeIn plus multiply object by a scalar
#' @export
setGeneric('mult', function(x, y, ...) standardGeneric('mult') )
#'
#' @describeIn mult multiply line by a scalar
#' @export
setMethod('mult', signature('line','numeric'),
          function(x, y, ...) {
            x <- l(x)
            l(cbind(x[,1:2,drop=FALSE]/y,x[,3]))
          } 
)
#' @details
#' matrix operates on the left on points in line as a column vector
#' @describeIn mult matrix operates on line as if multiplying the points in the line on the left 
#' @export
setMethod('mult', signature('line','matrix'),
          function(x, y, ...) {
            # matrix operates on the left on points in line as a column vector 
            x <- l(x)
            A <- solve(y)
            if(nrow(y) == 2 && ncol(y) ==2) A <- cbind(rbind(A, 0), c(0,0,1))
            l(x %*% A)
          } 
)
#' @export
setMethod('mult', signature('matrix','point'),
          function(x, y, ...) {
            # matrix operates on the left on points in line as a column vector 
            y <- p(y)
            A <- t(x)
            if(nrow(A) == 2 && ncol(A) ==2) A <- cbind(rbind(A, 0), c(0,0,1))
            # disp(y)
            # disp(A)
            # disp(y%*%A)
            p(y %*% A)
          } 
)
#' Rotation matrix
#' 
#' @param theta angle in radians that matrix would rotate
#' @examples
#'  library(spida2)
#'  rotm(pi/100)
#'  init()
#'  p(0,0,1) %>% pl
#'  l(1,1,1) %>% pl
#'  l(1,1,1) %>% mult(rotm(.1)) %>% pl
#'  l(1,1,1) %>% mult(rotm(pi/3)) %>% pl
#'  p(2,1,1) %>% {mult(rotm(pi/3),.)} %>% pl
#'  p(2,1,1) %>% {mult(rotm(0/3),.)} %>% pl
#'  
#'  
#'  l(1,1,0) %>% pl
#'  l(1,1,0) %>% mult(cbind(c(1,0),c(2,1))) %>% pl
#'  l(1,1,0) %>% mult(cbind(c(1,0),c(3,1))) %>% pl
#'  p(0,0,1) %>% pl(pch = 16)
#'  l(1,1,1) %>% pl
#'  l(1,1,1) %>% mult(cbind(c(1,0),c(3,1))) %>% pl
#'  l(1,1,1) %>% mult(cbind(c(0,-1),c(1,0))) %>% pl
#' @export
rotm <- function(theta = pi/2) {
  cbind( c(cos(theta), sin(theta)), c(-sin(theta), cos(theta)))
}
#' @export
setMethod('mult', signature('numeric','line'),
          function(x, y, ...) {
            mult(y, x, ...)
          } 
)

if(F){
  
  line2 %>% mult(2) %>% pl(col = 'red')
  mult(2, line2) %>% pl(col='blue')
  mult(-3, line2) %>% pl(col='blue')
  
  l(c(1,1,1))         
  l(c(1,1,1)) %>% mult(2) %>% pl  
}

#' @export
setMethod('mult', signature('point','numeric'),
          function(x, y, ...) {
            x <- p(x)
            p(cbind(x[,1:2, drop=FALSE]*y,x[,3]))
          } 
)
#' @export
setMethod('mult', signature('numeric','point'),
          function(x, y, ...) {
            mult(y, x, ...)
          } 
)

#' @export
setMethod('plus', signature('line','point'),
          function(x, y, ...) {
            obj <- l(x)
            pt <- p(y)
            
            nr <- nrow(obj)
            ind <- seq_len(nrow(pt))
            ind <- rep(ind, length.out = nr)
            pt <- pt[ind,,drop=FALSE]
            ret <- obj
            ret[,3] <- obj[,3,drop=F] - pt[,1,drop=F] * obj[,1,drop=F] - 
              pt[,2,drop =FALSE] * obj[,2,drop=FALSE]
            l(ret)
          } 
)

#' @export
setMethod('plus', signature('point','line'),
          function(x, y, ...) {
            plus(y, x, ...)
          } 
)

#' @export
setMethod('plus', signature('point','point'),
          function(x, y, ...) {
            obj <- p(x)
            pt <- p(y)
            nr <- nrow(obj)
            ind <- seq_len(nrow(pt))
            ind <- rep(ind, length.out = nr)
            pt <- pt[ind,,drop=FALSE]
            ret <- obj * pt[,3] + pt * obj[,3]
            ret[,3] <- ret[,3]/2
            p(ret)
          } 
)

#' @export
setMethod('plus', signature('gell','point'),
          function(x, y, ...) {
            x$center <- plus(p(x$center), p(y))[,1:2]
            x
          } 
)
#' @export
setMethod('plus', signature('point','gell'),
          function(x, y, ...) {
            plus(y, x, ...)
          } 
)

#' Drop a perpendicular
#' 
#' @param x,... line and point or point and line to drop perpendicular
#'              from point to line
#' @seealso [fta()]
#' @examples
#' init()
#' join(p(-1,-1),p(-2,2)) %>% pl -> line1
#' p(p(-1,-1),p(-2,2)) %>% pl
#' perp(line1) %>% pl(col = 'red')
#' p(0,0) %>% pl(pch=16)
#' perp(line1,p(2,2)) %>% pl(col = 'red')
#' 
#' @export
perp <- function(x,...) {
  UseMethod('perp')
}
#' @export
perp.line <- function(line, point = p(c(0,0))) {
  # line perpendicular through a given point
  line <- l(line)
  point <- p(point)
  ret <- line[,c(2,1,3), drop = FALSE]
  ret[,1] <- -ret[,1]
  ret[,3] <- 0
  plus(l(ret), point)
} 
#' @export
perp.point <- function(point, line) {
  line <-l(line)
  perp(line,point)
}
#' Turn each row of a matrix into one element of a list
#'
#' @param m matrix
#' @export
rows <- function(m) {
  # turn each row into element of list
  ret <- lapply(seq_len(nrow(m)), function(i) {
    ret <- m[i,,drop=FALSE]
    class(ret) <- class(m)
    ret
  })
  ret
}
#' Dual
#' 
#' @param x object of class ell, gell or ... to compute dual
#'
#' @export
dual <- function(x,...) {
  UseMethod('dual')
}
#' join now does the following:
# dual.line <- function(o1,o2) {
#   # line through 2 points or point at intersection of two lines
#   # - pipe through p or l if class uncertain
#   a1 <- rbind(o1)
#   a2 <- rbind(o2)
#   # recycle o2 if fewer rows
#   n1 <-nrow(a1)
#   n2 <-nrow(a2)
#   ind2 <- rep(seq_len(n2),length.out=n1)
#   ret <- pracma::cross(a1,a2[ind2,])
#   ret <- n(ret)
#   if('point' %in% class(o1)) class(ret) <- 'line'
#   if('line' %in% class(o1)) class(ret) <- 'point'
#   # if(!missing(class.out)) class(ret) <- class.out
#   ret
# }
#' @export
dual.ell <- function(x, ...) {
  shape <- attr(x,'params')$shape
  radius <- attr(x,'params')$radius
  ell(shape = solve(shape), radius = 1/radius)
}

#' @export
dual.gell <- function(x, ...) {
  gamma <- x$gamma
  d <- x$d
  radius <- x$radius
  gell(gamma = gamma, d= 1/d, radius = 1/radius)
}

if(FALSE){
  
  ge <- gell(A=cbind(c(2,1),c(0,3)))
  ge %>% pl
  ge %>% dual %>% pl
  gell() %>% pl(col = 'blue')
}
#' Generate points around a unit circle
#' 
#' @param n number of points

circle <- function(n = 100) {
  theta <- seq(0, 2*pi, length.out = n+1)
  cbind(cos(theta), sin(theta))
}

if(FALSE){
  init()
  circle() %>% cbind(1) %>% l  %>% pl
  circle() %>% p %>% mult(2) %>% pl(lwd = 3, type = 'l', col = 'red')
}

##' ell ####
#' probably no longer needed
# ell <- function (center = rep(0, 2), shape = diag(2), radius = 1, n = 100,
#                     angles = (0:n) * 2 * pi/n, fac = chol)
# {
#   rbindna <- function(x, ...) {
#     if (nargs() == 0)
#       return(NULL)
#     if (nargs() == 1)
#       return(x)
#     rbind(x, NA, rbindna(...))
#   }
#   sx <- c(shape[1,1] %>% sqrt, 0)
#   # disp(sx)
#   sy <- c(0, shape[2,2] %>% sqrt)
#   # disp(sy)
#   yx <- c(sx[1], shape[1,2] / sqrt(shape[1,1]))
#   # disp(yx)
#   xy <- c(shape[1,2] / sqrt(shape[2,2]), sy[2])
#   ry <- c(0, (shape[2,2] - shape[1,2]^2/shape[1,1]) %>% sqrt)
#   rx <- c((shape[1,1] - shape[1,2]^2/shape[2,2]) %>% sqrt, 0)
#   co <- c(sx[1],sy[2])
#   
#   # if (missing(ellipse) && missing(diameters) && missing(box))
#   #   all <- TRUE
#   circle <- function(angle) cbind(cos(angle), sin(angle))
#   rect <- cbind( c(-1,-1),c(-1,1),c(1,1), c(1,-1), c(-1,-1)) * sqrt(diag(shape))
#   Tr <- fac(shape)
#   ret <- t(c(center) + t(radius * circle(angles) %*% Tr))
#   params <- list(
#               parallelogram = t(c(center) + t(radius * rbind(c(1, 1), c(-1, 1), c(-1, -1),
#                                                              c(1,-1), c(1, 1)) %*% Tr)),
#               box = t(center + (radius * rect)),
#               yx = t(center + radius * cbind(yx, -yx)),
#               xy = t(center + radius * cbind(xy, -xy)),
#               diag1 = t(center + radius * cbind(co, -co)),
#               diag2 = t(center + radius * cbind(co*c(1,-1), -co*c(1,-1))),
#               ry = t(center + radius * cbind(ry, -ry)),
#               rx = t(center + radius * cbind(rx, -rx))
#   )
#   attr(ret,'params') <- c(params, list(center = center, shape = shape, radius = radius))
#   class(ret) <- 'ell'
#   ret
# }
# 
# params <- function(obj) {    # not needed ??????
#   attr(obj, 'params')
# }

## rot ####
#'
#' Rotate an object
#' 
#' @param m object to be rotated
#' @param theta angle by which to rotate object counterclockwise
#' 
#' @export 
rot <- function(m,theta,...) {
  UseMethod('rot')
}
#' @export 
rot.default <- function(m,theta) {
  ret <- m %*% cbind(c(cos(theta), -sin(theta)), c(sin(theta),cos(theta)))
  attributes(ret) <- attributes(m)
  ret
}
# rot.ell <- function(m,theta) {
#   params <- attr(m,'params')
#   rotm <- cbind(c(cos(theta), sin(theta)), c(-sin(theta),cos(theta)))
#   shape <- rotm%*%params$shape%*%t(rotm)
#   center <- params$center
#   radius <- params$radius
#   ell(center = center, shape = shape, radius = radius)
# }

#' @export 
rot.point <- function(m,theta) {
  ret <- m %*% cbind(c(cos(theta), -sin(theta),0), c(sin(theta),cos(theta),0),c(0,0,1))
  attributes(ret) <- attributes(m)
  ret
}
#' @export 
rot.line <- function(m,theta) {
  ret <- m %*% cbind(c(cos(theta), -sin(theta),0), c(sin(theta),cos(theta),0),c(0,0,1))
  attributes(ret) <- attributes(m)
  ret
}
#' @export 
rot.seg <- function(m,theta) {
  ret <- m %*% cbind(c(cos(theta), -sin(theta),0), c(sin(theta),cos(theta),0),c(0,0,1))
  attributes(ret) <- attributes(m)
  ret
}
#' Find intersection of ellipse with line
#'
#' @param ell an object of class ell
#' @param line an object of class line
#' 
#' @export
inter_ell_line <- function(ell, line) {
  # return 2 2-points
  center <- attr(ell,'params')$center %>% to2
  shape <- attr(ell,'params')$shape
  radius <- attr(ell,'params')$radius
  alphas <- chol(shape)%*%c(line)[1:2]
  theta <- atan2(-c(alphas)[1],c(alphas)[2])
  ret <- radius *t(chol(shape))%*%c(cos(theta),sin(theta))
  ret <- rbind(c(ret), -c(ret))
  ret
}
#' Find intersection of a generalized ellipse with line
#'
#' @param ge an object of class gell
#' @param line an object of class line
#' 
#' @export
inter_gell_line <- function(ge, line) {
  # return 2 2-points
  center <- ge$center
  gamma <- ge$gamma
  d <- ge$d
  radius <- ge$radius
  cline <- line %>% plus(p(center) %>% mult(-1))
  
  #########  need matrix times lines
  
  shape <- attr(ell,'params')$shape
  radius <- attr(ell,'params')$radius
  alphas <- chol(shape)%*%c(line)[1:2]
  theta <- atan2(-c(alphas)[1],c(alphas)[2])
  ret <- radius *t(chol(shape))%*%c(cos(theta),sin(theta))
  ret <- rbind(c(ret), -c(ret))
  ret
}


#'
#' ![ellipse intersections](Notes/Ellipse_intersections_centered.pdf)
#' 
# inter_ell <- function(ell1,ell2) {
#   # intersection of two centered ellipses
#   par1 <- attr(ell1,'params')
#   par2 <- attr(ell2,'params')
#   
#   cen <- sum(c(par1$center^2 ,par2$center^2))
#   if(cen > 0) warning('intersections are for ellipses centered at origin')
#   A <- par1$radius^2 * par1$shape
#   B <- par2$radius^2 * par2$shape
#   ell(shape = A) %>% pl
#   ell(shape = B) %>% pl(col='red')
#   
#   
#   T1 <- solve(t(chol(A)))
#   Bt <- T1%*%B%*%t(T1)
# 
#   eg <- eigen(Bt)
#   T2 <- t(eg$vectors) %*% T1
#   ell(shape = T2%*%A%*%t(T2)) %>% pl(lty=2)
#   ell(shape = T2%*%B%*%t(T2)) %>% pl(col='red')
#   lam1 <- eg$values[1]
#   lam2 <- eg$values[2]
#   disp(eg)
#   if(lam1 < 1 || lam2 > 1) {
#     return(c(NA,NA))
#   }
#   if(lam1 == lam2) {
#     return(c(NA,NA))
#   }
#   y2 <- lam2*(lam1 - 1)/(lam1 - lam2) 
#   x2 <- 1 - y2
#   ret <- rbind(
#     c(sqrt(x2),sqrt(y2)),
#     c(-sqrt(x2),sqrt(y2)),
#     c(-sqrt(x2),-sqrt(y2)),
#     c(sqrt(x2),-sqrt(y2))
#     
#   )
#   ret %>% points(pch=21:24, col = 'blue')
#   disp(ret)
#   disp(T1)
#   disp(T2)
#   ret <- ret %*% t(solve(T2))
#   ret
# }

# if(FALSE) {
#   Var <- .1*diag(2) + 1
#   ell1 <- ell(shape = Var)
#   ell2 <- ell(shape = diag(diag(Var)))
#   init()
#   inter_ell(ell1,ell2)  %>% points(pch = 21:24)
# 
#   ell1 %>% attr('params')$shape
#   
#     ell1 %>% pl(col='red')
#   ell2 %>% pl
#   
# }
#'
#' ## Generalized ellipses
#'
#' 
#' Constructor for gell objects
#' 
#' @param x an object of class ell to be transformed to a gell, or
#' @param center the center of the gell, default c(0,0), plus some combination of the following:
#' @param shape a \eqn{2 \times 2} matrix giving the variance corresponding to the ellipse
#' @param A a \eqn{2 \times 2} matrix with the transformation of the unit circle generating the shape
#'        of the ellipse
#' @param form a \eqn{2 \times 2} matrix with the matrix of the quadratic form, i.e. the inverse of `shape`
#' @param gamma an orthonormal matrix corresponding to the eigenvectors of 'shape'
#' @param d the eigenvalues of 'shape'
#' @param radius a factor by which to scale the ellipse as defined by other parameters 
#' @export
gell <- function(x,...) {
  if(missing(x)) gell.default(...)
  else UseMethod('gell')
}
setOldClass('gell')
#' @export
gell.ell <- function(ell,...) {
  params <- attr(ell,'parms')
  ei <- eigen(params$shape)
  ret <- list(gamma = ei$vectors, d = sqrt(ei$values), radius = params$radius, center = params$center)
  class(ret) <- 'gell'
  ret
}
setOldClass('ell')
#' @export
gell.default <- function(
    center = c(0,0), shape = var, var = solve(form), radius = 1, form = diag(2),
    A = NULL, gamma = eigen(shape)$vectors, d = sqrt(eigen(shape)$values)) {
  if(is.null(A)) {
    ret <- list(gamma = gamma, d = d, radius = radius, center = center)
  } else {
    sv <- svd(A, nv = 0)
    ret <- list(gamma = sv$u, d = sv$d, radius = radius, center = center)
  }
  class(ret) <- 'gell'
  ret
}

if(F){
  init()
  rbind(
    c(1,0,1),
    c(-1,0,1),
    c(0,-1,1),
    c(0,1,1)
  ) %>% l %>% pl
  gell() %>% pl
  ge <- gell(var = cbind(c(1,.6),c(.6,1)))
  ge %>% pl
  ge %>% dual %>%  pl
}


if(FALSE) {
  init()  
  gell(A=cbind(c(2,1),c(0,4))) %>%  pl(col = 'red')
  gell(A=cbind(c(2,1),c(0,4))) %>% rot(pi/4) %>% pl(col = 'red')
}

#' @describeIn gell the data ellipse of a data set with two variables contained in the 2xn matrix X
#' @export
gell.matrix <- function(X,...) {
  # data ellipse
  center <- apply(as.matrix(X), 2, mean)
  X <- t(t(X) - center)/sqrt(nrow(X))
  sv <- svd(X, nu = 0)
  ret <- list(gamma = sv$v, d = sv$d, radius = 1, center = center )
  class(ret) <- 'gell'
  ret
}  
#'
#' @describeIn pl method for objects of class 'ell'
#' @export
pl.ell <- function(ell,...) {
  pl(gell(ell), ...)
}
#' @describeIn pl method for objects of class 'ell'
#' @export
pl.gell <- function(gell, ...) {
  if(max(gell$d) < Inf & min(gell$d) >= 0) {
    U <- circle()
    pts <- t(gell$gamma %*% (t(U) * gell$radius * gell$d ) + c(gell$center))
    pl(seg(pts), ...)
  }
  else {
    # warning('pl.gell: case not yet implemented')
    fin <- gell$d != Inf
    fin1 <- gell$d[fin]*gell$radius*gell$gamma[,fin]+c(gell$center)
    fin2 <- -gell$d[fin]*gell$radius*gell$gamma[,fin]+c(gell$center)
    dir1 <- gell$gamma[,!fin]
    pts <- rbind(
      c(fin1,1),
      c(dir1,0),
      c(fin2,1),
      c(-dir1,0),
      c(fin1,1)
    )
    pts <- seg(pts)
    pl(pts,...)
  }
  invisible(gell)
}
#' Center of a gell object
#'
#' @export
cen <- function(x,...) {
  UseMethod('cen')
} 
#' center gell
#' 
#' @param gell a gell object
#' 
#' @export
cen.gell <- function(gell,...) {
  # generic in spida2
  p(c(gell$center))
}
#' Center of a gell object
#'
#' @export
center <- function(x,...) {
  UseMethod('center')
} 
#' center gell
#' 
#' @param gell a gell object
#' 
#' @export
center.gell <- function(gell,...) {
  # generic in spida2
  p((gell$center))
}
#' Next
#' 
#' @importFrom spida2 center
#' @export
pcenter <- center

if(FALSE) {
  
  X <- cbind(rnorm(200),rnorm(200)) %*% cbind(c(1,1), c(.2,-.3))
  init()
  X %>% p %>% pl  
  gell(X) %>% pl
  methods(pl)
  
  gell(X) %>% center.gell %>% pl(pch = 16, col = 'red')
  
}



## axes ####
#' S4 method to find axes of a generalized ellipse
#' 
#' @export
setGeneric('axes', function(x,y,type,tan,end,...) standardGeneric('axes'))
setOldClass('gell')
setOldClass('point')

## axes: gell point missing missing missing ####

#' @export
setMethod("axes", signature(x = 'gell', y = 'point', type = 'missing', tan = 'missing', end = 'missing'),
          function(x, y, type, ...)  {
            P <- cbind(c(0,1),c(-1,0))
            pt <- p(y)
            gell <- x
            pt <- plus(pt, mult(pcenter(gell),-1)) # centered ellipse
            v <- diag(1/gell$d) %*% t(gell$gamma) %*% cbind(pt[1:2])
            u <- v / sqrt(sum(v^2))
            pts <- gell$radius * gell$gamma %*% diag(gell$d) %*% cbind(u, P%*%u)
            pts_on_ell <- p(t(pts)) %>% plus(center(gell))
            tan_dirs <- p(t(pts))[2:1,]
            tan_dirs[,3] <- 0
            list(points = pts_on_ell, tan_dirs = p(tan_dirs))
          })



## axes: gell point character logical missing ####

#' @export
setMethod("axes", signature(x = 'gell', y = 'point', type = 'character', tan = 'logical', end = 'missing'),
          function(x, y, type, tan, ...)  {
            # types:
            #     l - line
            #       tan F: point is direction
            #       tan T: point is tangent direction to ellipse
            #     h - half line from center
            #     r - radius
            #     d - diameter
            library(spida2)
            ax <- axes(x,y)
            if(type == 'l' ){  # line
              if(!tan) {
                ret <- join(center(x),p(ax$points[1,]))
              } else {
                ret <- join(center(x),p(ax$points[2,]))
              }
            }
            if(type == 'h') {  # half-line
              if(!tan) {
                ret <- rbind(
                  center(x),
                  ax$tan[2,]
                ) %>% seg
              } else {
                ret <- rbind(
                  center(x),
                  ax$tan[1,]
                ) %>% seg
              }
            }
            if(type == 'r') {   # radius
              if(!tan) {
                ret <- rbind(
                  center(x),
                  ax$points[1,]
                ) %>% seg
              } else {
                ret <- rbind(
                  center(x),
                  ax$points[2,]
                ) %>% seg
              }
            }
            if(type == 'd') {   # diameter
              if(!tan) {
                ret <- rbind(
                  plus(p(center(x)),plus(p(center(x)), mult(p(ax$points[1,]),-1))),
                  ax$points[1,]
                ) %>% seg
              } else {
                ret <- rbind(
                  plus(p(center(x)),plus(p(center(x)), mult(p(ax$points[2,]),-1))),
                  ax$points[2,]
                ) %>% seg
              }
            }
            ret
          }
)

## axes: gell point character logical line ####

#' @export
setMethod("axes", signature(x = 'gell', y = 'point', type = 'missing', tan = 'logical',end = 'line'),
          function(x, y, type, tan, end){
            li <- axes(x, y, type = 'l', tan = tan)
            ret <- rbind(
              center = center(x),
              end = join(li, end))
            seg(p(ret))
          }
)


## axes: gell point missing missing line ####

#' @export
setMethod("axes", signature(x = 'gell', y = 'point', type = 'missing', tan = 'missing',end = 'line'),
          function(x, y, type, tan, end){
            axes(x, y, tan = FALSE, end = end)
          }
)

## axes: gell point character missing missing ####
#' @export
setMethod("axes", signature(x = 'gell', y = 'point', type = 'character', tan = 'missing',end = 'missing'),
          function(x, y, type, tan){
            axes(x, y, type, FALSE)
          }
)



if(F){
  
  init()
  ge %>% pl
  axes(ge, p(1,0,0), 'l')
  
  ge %>% axes(p(1,0,0), type = 'l', tan = F) %>% pl %>% print
  ge %>% axes(p(1,0,0), end = l(1,0,0), tan = F) %>% pl(col=green,lwd =2) %>% print   
  ge %>% axes(p(1,0,0), type = 'l',F) %>% pl
  ge %>% axes(p(1,0,0), type = 'l',T) %>% pl
  ge %>% axes(p(1,0,0), type = 'r',F) %>% pl(col=red, lwd = 2)
  ge %>% axes(p(1,0,0), type = 'r',T) %>% pl(col=red, lwd = 2)
  ge %>% axes(p(1,0,0), type = 'd',F) %>% print %>% pl(col=blue, lwd = 2)
  ge %>% axes(p(1,0,0), type = 'd',T) %>% print %>% pl(col=blue, lwd = 2)
  ge %>% axes(p(1,0,0), type = 'h',F) %>% pl(col=green, lwd = 2)
  ge %>% axes(p(1,0,0), type = 'h',T) %>% pl(col=green, lwd = 2)
  
  for(i in 0:20){
    ge %>% axes(rot(p(1,0,0),i*pi/10), type = 'h',F) %>% pl(col=green, lwd = 2)
  }
}

#' @export
text.point <- function(pp, labels = seq_len(nrow(pp)), ...) {
  # adj:  x,y in [0,1]
  # pos: overrides adj: 1:4 : below,left,above,right
  # offset: with pos in character widths
  # cex, srt, family, xpd
  
  pp <- p(pp)
  ret <- pp
  pp[pp[,3] == 0,] <- NA 
  # labels <- TeX(labels)
  text(pp[,1], pp[,2], labels = labels, ...)
  invisible(ret)
}

#' @export
sel <- function(x, what, ...) {
  UseMethod('sel')
}

#' @export
sel.list <- function(x, what,...) {
  x[[what]]
}

#' @export
sel.default <- function(x, what,...) {
  if(length(dim(x)) == 2) ret <- x[what,,drop=FALSE]
  class(ret) <- class(x)
  ret
}

if(F) rbind(c(1,1,0), c(2,1,0)) %>% p %>% sel(1) 


setOldClass('gell')
setOldClass('point')
setOldClass('line')

if(FALSE) {
  library(latex2exp)
  X <- matrix(rnorm(20), ncol = 2)
  X <- t(c(1,2)+t(X))
  init()
  X %>% pl
  gell(X) %>% pl
  gell(X) %>% axes(p(1,1), type='l') %>% pl(col = 'red')
  gell(X) %>% axes(p(1,1)) %>% sel('tan_dirs') %>% join(center(gell(X))) %>% pl
  gell(X) %>% axes(p(1,1), end = l(0,1,0))  %>% pl ->li
  li[2,] %>% p %>% text(TeX('$\\beta_1$'), pos = 1)
  gell(X) %>% gell_box(p(1,0,0),p(0,1,0))
  pl(col = 'gray')
  gell(X) %>% center %>% pl(col = 'red')
  gell(X) %>% axes(p(1,0,0)) %>% pl
}

# gell_conj <- function(ge, dir) {
#   # line thru center conjugate to a direction
#   # i.e. intersecting ellipse with tangent parallel to a direction
#   a <- axes(ge, dir)
# }

#' @export
gell_box <- function(ge, p1, p2 = p(a1$points[2,])) {
  # works: %>% 2023-02-24
  a1 <- axes(ge, p1)
  l1 <- join(p(a1$points[2,,drop=FALSE]), p(a1$tan_dirs[2,,drop=FALSE]))
  rad <- p(a1$points[2,]) %>% plus(center(ge) %>% mult(-1))
  off <- rad %>% mult(-2)
  l3 <- plus(l1,off) 
  a2 <- axes(ge, p2)
  l2 <- join(p(a2$points[2,,drop=FALSE]), p(a2$tan_dirs[2,,drop=FALSE]))
  rad <- p(a2$points[2,]) %>% plus(center(ge) %>% mult(-1))
  off <- rad %>% mult(-2)
  l4 <- plus(l2,off) 
  l(rbind(l1,l2,l3,l4))
}

#' @export
gell_box_ <- function(ge) {
  gell_box(ge, p(0,1,0), p(1,0,0))
}

if(F){
  init()
  l(1,1,1) %>% pl
  l(1,1,1) %>% plus(p(1,1,1))%>% pl
  p(0,0,1) %>% pl(pch = 16)
  p(1,1,1) %>% pl(pch = 16)
  
  gell(X) %>% pl
  gell(X) %>% gell_box(p(1,0,0),p(0,1,0)) %>% pl
  gell(X) %>% gell_box(p(1,0,0),p(0,1,0)) %>% seg %>% close %>% pl(fill='#00009922')
  gell(X) %>% rot(pi/7) %>% seg %>% pl(fill= "#44000044", border = 'blue')
  gell(X) %>% gell_box(p(1,1,0)) %>% pl
  gell(X) %>% gell_box(p(1,1,0)) %>% seg %>% pl
  gell(X) %>% gell_box(p(1,0,0),p(0,1,0)) %>% pl
  gell(X) %>% gell_box(p(1,0,0)) %>% pl(col = 'red')
  gell(X) %>% gell_box(p(0,1,0)) %>% pl(col = 'green')
  
  
  gell(X) %>% center %>% pl(col = 'red') %>% {.[1,]} %>% p %>% pl
  gell(X) %>% axes(p(1,1,0)) %>% pl(col='red', pch = 16) %>% print
  join(p(1,1,0), center(gell(X))) %>% pl(col = 'red') 
  join(p(1,1,1), center(gell(X))) %>% pl(col = 'green')
  
  gell(X) %>% axes(p(1,1,0)) %>% pl(col='red', pch = 16) %>% print
  gell(X) %>% axes(p(1,1,1)) %>% pl(col='green', pch = 16) %>% print
  
  gell(X) %>% center %>% pl
  gell(X) %>% gell_box(p(1,1,0)) %>% pl
  gell(X) %>% axes(p(-c(0,1))) %>% pl(pch = 16, col = 'red')
  gell(X) %>% axes(-c(0,1)) %>% pl(pch = 16, col = 'red')  
}




# axes.gell <- function(gell, obj, ...) {
#   axes_gell(obj, gell, ...)
# }
# axes_gell <- function(obj, gell, ...) {
#   UseMethod('axes_gell')
# }
# axes_gell.default <- function(obj, gell, ...) {
#   pt <- p(obj)
#   pt <- plus(pt, mult(pcenter(gell),-1)) # centered ellipse
#   v <- diag(1/gell$d) %*% t(gell$gamma) %*% cbind(pt[1:2])
#   u <- v / sqrt(sum(v^2))
#   pt <- gell$gamma %*% diag(gell$d) %*% u
#   p(c(pt)) %>% plus(pcenter(gell))
# }


#######   Intersections     #####################


#' @export
inter_gell_line <- function(ge, line) {
  # GOOD: I THINK
  inter_circle_line <- function(line) {
    rotm <- function(theta = pi/2) {
      cbind( c(cos(theta), -sin(theta)), c(sin(theta), cos(theta)))
    }
    hypot <- pracma::hypot
    line <- l(line)
    # work out intersections if line has form (h, 0, 1) or (h, 0, 0) i.e. y-axis
    h <- hypot(line[,1]/line[,3],line[,2]/line[,3])
    # If h < 1, no intersection, if h is NaN, line is y axis
    x <- - ifelse(is.na(h), 0, 1/h) # is.na(x) implies the y-axis
    y <- sqrt(1 - x^2) # is NaN if no intersection
    p1 <- cbind(x,y,1)
    p1[is.na(y),] <- 0
    p2 <- p1
    p2[,2] <- -p2[,2]
    theta <- atan2(line[,2],line[,1])
    rots <- lapply(theta, rotm)
    for(i in seq_along(rots)) {
      p1[i,1:2] <- p1[i,1:2,drop=FALSE] %*% rots[[i]] 
      p2[i,1:2] <- p2[i,1:2,drop=FALSE] %*% rots[[i]] 
    }
    list(p(p1),p(p2))
  }
  
  line <- line %>% plus(mult(-1,center(ge)))
  line <- line %>% mult(diag(1/ge$d)%*%t(ge$gamma))
  ps <- inter_circle_line(line)
  ps <- lapply(ps, function(pt) {
    mult(ge$gamma %*% diag(ge$d), pt) %>% plus(center(ge))
  })
  ps
} 
if(F){
  # play with axes to solve intesection problem
  ge <- gell(X)
  init()
  ge %>% pl %>% center %>% pl
  join(center(ge),p(1,0,0)) %>% pl
  join(center(ge),p(-1,0,0)) %>% pl(col= red)
  axes(ge,p(1,0,0)) %>% pl
  axes(ge,p(-1,0,0)) %>% pl(col = blue)
  axes(ge,p(1,1,0)) %>% pl(col = green)
  
  # seems: axes gives the right point in direction of point at infty
}

#' @export
inter_gell_seg <- function(ge, sg) {
  # %>% HIGHLY %>% BROKEN
  # maybe we won't need this is including circle includes all finite points
  # because we can then use axes to find intersection from interior to infinite point
  if(nrow(sg) == 1) sg <- seg(rbind(center(ge), sg))
  if(nrow(sg) >2) warning('only first two rows of sg are used')
  li <- join(p(sg[1,]), p(sg[2,])) # assuming distinct
  ps <- inter_gell_line(ge, li)        # intersection with line joining two points
  ps <- p(rbind(sg,ps[[1]], ps[[2]]))  # combine
  ps <- plus(ps, mult(p(ps[3,]), -1))  # center on one   to find order
  
  pt2 <- ps[2,1:2]                       
  mat <- rbind( c(pt2), c(-pt2[2],pt2[1]))/(sum(pt2^2))
  mat <- rbind( cbind(mat,0), c(0,0,1))
  ps <- mult(mat, ps)  # TEST on set including pt at infty ##########################
  
  line <- line %>% plus(mult(-1,center(ge)))
  line <- line %>% mult(diag(1/ge$d)%*%t(ge$gamma))
  ps <- inter_circle_line(line)
  ps <- lapply(ps, function(pt) {
    mult(ge$gamma %*% diag(ge$d), pt) %>% plus(center(ge))
  })
  ps
} 

#' 
#'  From @alberich-carraminana2017
#' 

#' @export
disp <- spida2::disp
#'
#' Plot contour lines for a function
#' 
#' Uses [graphics::contour] to draw contour lines
#' 
#' @param function of two arguments
#' @param val numerical values defining contour lines of functions, default: NULL
#' @param xlim 
#' @param ylim 
#' @param fac factor by which to expand `par("usr")` to create grid
#' @param n resolution of grid to define contours, default: 200
#' @param ... arguments passed to plotting function, e.g. col, lty, ... 
#' @returns adds contour lines to existing plot
#' 
#' @export 
contour_function <- function(fun, val = NULL, xlim = expand(par("usr")[1:2], fac), 
                             ylim = expand(par("usr")[3:4], fac), fac = 1.1, n = 200, ...) {
  cont <- function(x, y, fun, ..., nlevels = 50, levels = NULL) {
    ret <- as.matrix(expand.grid(x, y))
    ret <- apply(ret, 1, function(v) fun(v[1], v[2]))
    ret <- matrix(ret, nrow = length(x), ncol = length(y))
    if (is.null(levels)) {
      graphics:::contour.default(x, y, ret, nlevels = nlevels, add = TRUE, ...)
    }
    else {
      graphics:::contour.default(x, y, ret, levels = levels, add = TRUE, ...)
    }
  }
  expand <- function(v, fac) {
    cen <- mean(v)
    cen + fac * (v - cen)
  }
  x <- seq(xlim[1], xlim[2], length.out = n)
  y <- seq(ylim[1], ylim[2], length.out = n)
  cont(x, y, fun, levels = val, ...)
}
#' @describeIn contour_function method for functions for contour generic function
#' @export
contour.function <- contour_function
#' @export
inter_fun <- function(fun, val = 0, 
                      xlim = expand(par('usr')[1:2],fac),
                      ylim = expand(par('usr')[3:4],fac),
                      fac = 1.1, n = 200) {
  cont <- function(x,y,fun, ..., nlevels = 50, levels=NULL) {
    ret <- as.matrix(expand.grid(x,y)) 
    # print(head(ret))
    ret <- apply(ret, 1, function(v) fun(v[1],v[2],...))
    ret <- matrix(ret, nrow = length(x), ncol = length(y))
    # print(dim(ret))
    if(is.null(levels)) {
      ret <- contourLines(x,y,ret, nlevels = nlevels)
      
    } else {
      ret <- contourLines(x,y,ret, levels = levels)
    }
  }
  expand <- function(v, fac) {
    cen <- mean(v)
    cen + fac*(v - cen)
  }
  x <- seq(xlim[1], xlim[2], length.out = n)
  y <- seq(ylim[1], ylim[2], length.out = n)
  conts <- cont(x, y, fun, levels = val)
  conts <- lapply(conts,
                  function(x) {
                    list(level = x$level,
                         contour = seg(cbind(x$x,x$y)))
                  })
  class(conts) <- 'contour_list'
  conts
}
if(F){
  
init()
fun <- function(x,y) x^3 - y^2
inter_fun(fun,val = c(-1,-.1,0,.1,1), n =400) %>% pl()
}

#' @export
pl.contour_list <- function(cl,...) {
  lapply(cl, function(x) {
    pl(x$contour, ...)
  }
  )
  invisible(cl)
}
#' Intersecting ellipses
#' 
#' Normally called from S4 method [inter]
#' 
#' @param g1 gell object
#' @param g2 gell object
#' @param n fineness of grid to compute intersections
#' @param verbose default: 0
#' @export
inter_gells <- function(g1, g2, n = 1000, verbose = 0) {
  # intersection points of two ellipses
  if(verbose>0 || options('verbose')$verbose) {
    init()
    abline(v=0)
    abline(h=0)
    g1 %>% pl(col='red')
    g2 %>% pl(col='red', lty = 2)
  }
  # 
  # transform shape of g1 to unit circle
  # 
  T2C <- (1/g1$radius) * diag(1/g1$d) %*% t(g1$gamma)
  #
  # Apply transform to g2
  #
  g2c <- gell(gamma = g2$gamma, radius = g2$radius, d = g2$d, center = g2$center - g1$center)
  gt2 <- gell(A = T2C %*% g2c$gamma %*% diag(g2c$d), center = T2C %*% g2c$center, radius = g2c$radius)
  disp(gt2)
  if(verbose > 0) {
    gt2 %>% pl(col = 'green',lty = 2)
    gell() %>% pl(col = 'green')
  }
  #
  # Matrix to transform unit circle to shape of gt2
  #
  A2 <- gt2$radius * gt2$gamma %*% diag(gt2$d)
  #
  # Function whose zeros identify points in gt2 that are on unit circle
  #
  fun <- function(theta) {
    sum((A2 %*% rbind(sin(theta), cos(theta)) + gt2$center)^2) - 1
  }
  #
  # Vectorized version for uniroot.all
  # 
  funv <- function(theta) {
    sapply(theta, fun)
  }
  roots <- rootSolve::uniroot.all(funv, c(0,2*pi), n = n )
  if(verbose>0)disp(roots)
  if(length(roots) > 0){
    points <- lapply(roots,
                     function(rt)
                       A2 %*% rbind(sin(rt), cos(rt)) + gt2$center
    )
    if(verbose>0)disp('green')
    if(verbose>0)disp(points)
    if(verbose > 0) p(t(do.call(cbind,points))) %>% pl(col='green')
    # transform to points on g1/g2
    points <- lapply(points,
                     function(p) {
                       solve(T2C) %*% p + g1$center
                     })
    ret <- p(t(do.call(cbind,points)))
    if(verbose >0) pl(ret, pch = 16, col = 'red')
  } else {
    ret <- numeric(0)
    if(verbose > 0) cat("\nno intersection found\n")
  }
  ret
}
#' Generate the intersection of two objects
#' 
#' @param x gell or function
#' @param y gell, line, segment, point or numeric
#' @param option fineness of grid over which to find intersections, default: 
#' 
#' @returns intersecting point(s)
#' @export
setGeneric('inter', function(x, y, option,...) standardGeneric('inter'))

## inter: gell gell missing ####
#' @describeIn inter intersection of two `gell` objects using [inter_gells]
#' @export
setMethod("inter", signature(x = 'gell', y = 'gell', option = 'missing'),
          function(x, y, option, ...)  {
            inter_gells(x, y)
          }
)

## inter: gell gell numeric ####
#' @describeIn inter intersection of two `gell` objects
#' @export
setMethod("inter", signature(x = 'gell', y = 'gell', option = 'numeric'),
          function(x, y, option, ...)  {
            inter_gells(x, y, n=option)
          }
)

## inter: gell line missing ####
#' @describeIn inter intersection a `gell` object with a `line`
#' @export
setMethod("inter", signature(x = 'gell', y = 'line', option = 'missing'),
          function(x, y, option, ...)  {
            inter_gell_line(x, y)
          }
)

## inter: gell seg missing ####
#' @describeIn inter intersection a `gell` object with a `seg` object (segment)
#' @export
setMethod("inter", signature(x = 'gell', y = 'seg', option = 'missing'),
          function(x, y, option, ...)  {
            inter_gell_seg(x, y)
          }
)

## inter: gell point missing ####
#' @describeIn inter intersection a `gell` object with a `point`
#' @export
setMethod("inter", signature(x = 'gell', y = 'point', option = 'missing'),
          function(x, y, option, ...)  {
            inter_gell_seg(x, p(y))
          }
)

## inter: function numeric ####
#' @describeIn inter intersection a `function` object with a numeric value(s)
#' @export
setMethod("inter", signature(x = 'function', y = 'numeric', option = 'missing'),
          function(x, y, option, ...)  {
            inter_fun(x, y,...)
          }
)


#' @export
rmat <- function(rho) {
  cbind(c(1,rho), c(rho, 1))
}

if(F){
  
  g1 <- gell(shape=rmat(-.8), center = c(-1,1))
  g2 <- gell(shape=rmat(.9), radius = 2, center = c(-1,1))
  gs <- gell(shape=cbind(c(1,-1), c(-1,1)), radius = 2, center = c(-1,1))
  
  
  init()
  gs %>% pl
  
  inter(g1,g2) %>% pl
  init()
  g1 %>% pl
  g2 %>% pl
  
  gell(shape=cbind(c(1,0),c(0,1))) %>% pl -> gc
  gell(shape=cbind(c(1,0),c(0,.5)), radius = 1) %>% pl -> gf
  gell(A=rotm(pi/6)%*%cbind(c(1,0),c(0,1))) %>% pl -> gc
  gell(A=rotm(pi/6)%*%cbind(c(1,0),c(0,.5)), radius = 1.00000000001) %>% pl -> gf
  inter_gells(gc,gf) %>% pl %>% print
  # debug(inter_gells)
  
  
  
  X %>% pl
  init()
  gell(X) %>% pl
  gell() %>% pl
  l(-.5,-.5,1) %>% pl -> l1
  l(-.2,-.5,1) %>% pl -> l2
  l(-.5,-.5,0) %>% pl -> l3
  ls <- rbind(l1,l2,l3) %>% l
  
  gell(X) %>% inter_gell_line(ls) %>% pl
  
  init()
  gell() %>% pl -> ge
  p(0,0,1) %>% pl
  
  
  
  l(1,1,1) %>% pl
  
  
  # let's see whether points at infinity work as they should
  # 
  X <- t(c(2,3) + t(X))
  
  X %>% plot
  gell(X) %>% pl(col = 'red')
  gell(X) %>% axes(p(c(1,0))) %>% pl(pch = 16, col = 'red')
  gell(X) %>% pcenter %>% pl(pch = 16, col = 'blue')
  gell(X) %>% axes(p(1,0,0)) %>% pl(pch = 16, col = 'red')
  gell_box(gell(X), p(1,0,0)) %>% pl
  gell_box(gell(X), p(1,0,0), p(0,1,0)) %>% pl
}

#' @export
dual.gell <- function(gell, center = gell$center) {
  ret <- gell
  ret$d <- 1/ret$d
  ret$center <- center
  ret
}

if(F){
  
  X %>% plot(asp = 1)
  gell(X) %>% pl(col = 'red')
  
  gell(X) %>% dual %>% pl
  (center = pcenter(gell(X))[[1:2]])
  
}

if(F){
  init()
  rbind(
    c(1,1,1),
    c(1,2,2),
    c(1,-1,0),
    c(2,2,1)
  ) %>% p %>% seg -> s1
  s1
  init()
  s1 %>% pl(close=T, fill = red) -> so1
  att <- attr(so1,'att')  
  att %>% lapply(class)
  att %>% lapply(pl)
  
  init()
  rbind(
    c(1,1,1),
    c(1,2,1),
    c(-2,-2,1)
  )  %>% seg %>% close %>% pl(fill = "#0000FF22")
  init()
  rbind(
    c(1,-1,0),
    c(1,2,1),
    c(1,0,0),
    c(-2,-2,1),
    c(-1,0,0)
  )  %>% seg %>% rot(pi/4) %>% close %>% pl(fill = "#0000FF22")
  
  init()
  ge <- gell(X)
  ge %>% pl
  ge %>% seg %>% pl(fill = "#0000FF22")
  gz <- gell(d = c(Inf,0))
  gz %>% seg %>% pl
  
  redish <- "#88000033"
    rbind(
      c(-1,-1,1),
      c(1,-1,0),
      c(1,1,1),
      c(-1,1,0),
      c(-1,-1,1)
    )  %>% seg %>% close %>% rot(2*pi/6) %>% pl(fill = redish)
    init()
    rbind(
      c(1,-1,0),
      c(-1,-1,1),
      c(-1,1,0),
      c(1,1,1),
      c(1,-1,0)
    )  %>% seg %>% close %>% rot(pi/16)%>% pl(fill = "#FF00FF22")
    
}
if(F){
  
  p(l(c(1,1,1)),l(c(2,1,1)))
  l(c(1,1,1)) %>% pl
  l(c(2,1,1)) %>% pl
  p(l(c(1,1,1)),l(c(2,1,1))) %>% pl(pch = 16)
  p(l(c(1,1,1),c(2,1,1))) %>% pl(col = 'red', lwd = 2)
  p(l(c(1,1,1))) %>% pl(col = 'red', lwd = 2)
  l(c(1,1,3),c(2,1,3)) %>% pl(col = 'green')
  
  init()
  l(c(1,1,3),c(2,1,3),c(1,1,0)) %>% pl(col = 'cyan')
  l(c(1,1,3),c(2,1,3)) %>% p(l(c(1,1,0))) %>% pl(col = 'purple')
  
  
}


if(F){
  
  init()
  zell %>% pl
  zell %>% dual %>% pl
  zell %>% rot(.1) %>% pl(col='red')
  lapply(seq(0,pi,.2), function(r) zell %>% rot(r) %>% pl(col='red')) %>% invisible
  zell %>% params 
  
  zl <- l(c(1,1,0))
  zl %>% pl
  zl %>% rot(.1) %>%  pl
  zp <- p(c(1,-1,1))
  zp %>% pl(pch=16,col='red')
  zp %>% rot(.1) %>% pl(pch=16,col='red')
  zp %>% plus(c(1,1)) %>% pl(pch=16,col='red')
  
  
  init()
  seg(rbind(c(1,1,0),c(1,2,1),c(3,3,1)) %>% l) %>% pl
  rbind(c(1,1,0),c(1,2,1),c(3,1,1)) %>% l %>% pl
  rbind(c(1,1,0),c(1,2,1),c(3,3,1)) %>% l %>% pl
  
  p(matrix(rnorm(30),ncol=3))
  pl(z)
  
  
  init()
  
  
  ell(shape = .1*diag(2) + 1, radius = 2) %>% pl(col='red')
  ell(shape = .1*diag(2) + 1, radius = 2) %>% dual %>% pl
  ell(shape = .1*diag(2) + 1) %>% dual %>% pl
  
  
  
}

## from to along ##############    

## fta ####

#' @export
setGeneric('fta', function(x,y,z,...) standardGeneric('fta'))

## fta point point missing ----------
#' @export
setMethod("fta", signature(x = 'point', y = 'point', z = 'missing'),
          function(x, y, z, ...)  {
            seg(rbind(x,y))
          }
)

## fta point line missing ----------
#' @export
setMethod("fta", signature(x = 'point', y = 'line', z = 'missing'),
          function(x, y, z, ...)  {
            to <- join(perp(y,x),y)
            seg(rbind(x,to))
          }
)

## fta point point missing ----------
#' @export
setMethod("fta", signature(x = 'line', y = 'point', z = 'missing'),
          function(x, y, z, ...)  {
            from <- join(perp(x,y),x)
            seg(rbind(from,y))
          }
)

## fta line line line ----------
#' @export
setMethod("fta", signature(x = 'line', y = 'line', z = 'line'),
          function(x, y, z, ...)  {
            # z is used as a direction
            from <- join(x,z)
            to <- join(y,z)
            seg(rbind(from,to))
          }
)

## fta line line point ------------------------
#' @export
setMethod("fta", signature(x = 'line', y = 'line', z = 'point'),
          function(x, y, z, ...)  {
            stop('no method from line to line along point')
          }
)

## fta line point line ------------------------------
#' @export
setMethod("fta", signature(x = 'line', y = 'point', z = 'line'),
          function(x, y, z, ...)  {
            # z is used as a direction
            z[3] <- 0
            z <- l(z)
            z <- plus(z,y)
            from <- join(x, z)
            to <- y
            seg(rbind(from,to))
          }
)

## fta point line line ---------------------------
#' @export
setMethod("fta", signature(x = 'point', y = 'line', z = 'line'),
          function(x, y, z, ...)  {
            ret <- fta(y, x, z)
            seg(ret[2:1,])
          }
)


## fta line point point ---------------------------
#' @export
setMethod("fta", signature(x = 'line', y = 'point', z = 'point'),
          function(x, y, z, ...)  {
            # z used as a direction
            z <- join(z,p(0,0,1))
            fta(x,y,z)
          }
)

## fta point line point  -------------------------------
#' @export
setMethod("fta", signature(x = 'point', y = 'line', z = 'point'),
          function(x, y, z, ...)  {
            z <- join(z,p(0,0,1))
            fta(x,y,z)
          }
)
if(F){
  
  init()
  gell(shape = rmat(.8), center = c(1,1)) %>% pl -> ge
  ge %>% axes(p(0,1,0)) %>% pl
  ge %>% axes(p(0,1,0)) %>% sel(2) %>%  join(p(0,0,1)) %>% pl
  ge %>% axes(p(0,1,0)) %>% {join(.$points,.$tan_dirs)} %>% pl
  ge %>% axes(p(0,1,0))
  
}


## fta gell line point -------------------------------------
#' @export
setMethod("fta", signature(x = 'gell', y = 'line', z = 'point'),
          # use point a direction towards line
          # so only the center matters
          function(x, y, z, ...)  {
            from <- p(x$center)
            fta(from,y,z)
          }
)

## fta gell line line -------------------------------------
#' @export
setMethod("fta", signature(x = 'gell', y = 'line', z = 'line'),
          # use z as tangent to direction
          # simply use l(x,y,0) to specify tangent direction
          function(x, y, z, ...)  {
            along <- p(z[2], -z[1], 0)
            along <- axes(x, along, 'l', tan=T)
            from <- p(x$center)
            fta(from,y,along)
          }
)

## fta gell point missing -------------------------------------
#' @export
setMethod("fta", signature(x = 'gell', y = 'point', z = 'missing'),
          function(x, y, z, ...)  {
            fta(p(x$center),y)
          }
)

#' @export
comb <-function(g1,g2,lam = .5) {
  w1 <- (1-lam)  * g1$gamma %*% diag(1/g1$d^2) %*% t(g1$gamma)/g1$radius
  w2 <- (lam) * g2$gamma %*% diag(1/g2$d^2) %*% t(g2$gamma)/g2$radius
  cen <- w1 %*% g1$center + w2 %*% g2$center
  var <- solve(w1 + w2)
  cen <- var %*% cen
  gell(shape = var, center = cen)
}

## fta gell gell missing
#' @export
setMethod("fta", signature(x = 'gell', y = 'gell', z = 'missing'),
          function(x, y, z, ...)  {
            ret <- sapply(seq(0,1, length.out = 1001) , 
                          function(lam) comb(x, y, lam)$center)
            seg(p(t(ret)))
          }
)
#' @export
setMethod("fta", signature(x = 'gell', y = 'gell', z = 'numeric'),
          function(x, y, z, ...)  {
            ret <- sapply(seq(0,1, length.out = z) , 
                          function(lam) comb(x, y, lam)$center)
            seg(p(t(ret)))
          }
)

if(F){
  init()
  fta(g1,g2,5)  %>% l %>% pl(pch = 16,cex = 2 , col = 'red')
  fta(center(g1),l(0,1,-4))  %>% l %>% pl(pch = 16,cex = 2 , col = 'red')
  l(0,-1,4) %>% fta(fta(g1,g2,5) %>% p,.) %>% pl
}


## ptr generic for projective transformations  -------------------------
#' @export
setGeneric('ptr', function(x,y,z,...) standardGeneric('ptr'))

## ptr point matrix missing ----------
#' @export
setMethod("ptr", signature(x = 'point', y = 'matrix', z = 'missing'),
          function(x, y, z, ...)  {
            p(t(y %*% t(x)))
          }
)

## ptr line matrix missing ----------
#' @export
setMethod("ptr", signature(x = 'line', y = 'matrix', z = 'missing'),
          function(x, y, z, ...)  {
            l(t(t(y) %*% t(x)))
          }
)

setOldClass('conic')

conic <- function(x,...) {
  useMethod('conic')
}

conic.gell <- function(x,...) {
  
}


## Examples ########################################

if(F){
  library(spida2)
  init()
  abline(h=0)
  abline(v=0)
  gell(shape = rmat(.8), c = c(1,2)) %>% pl -> ge
  l(1,1,-1) %>% pl
  l(-1,-2,1) %>% pl(col = 'cyan')
  p(-1,-2,1) %>% pl(col = 'cyan', pch = 16)
  p(-1,-2,1) %>% join(p(0,0,1))%>% pl(col = 'cyan', pch = 16)
  fta(ge, l(1,1,-1), p(-1,-2,1)) %>% pl
  p(1,1,0) %>% pl
  axes(ge, p(1,2,0), 'l', T) %>% pl
  
  fta(ge, l(1,1,-1), l(-1,-2,1)) %>% pl
  l(1,1,-1) %>% pl
  
  fta(ge, l(1,1,-1), l(0,1,0)) %>% pl
  fta(ge, p(0,0,1)) %>% pl
  fta(ge, p(1,-1,0)) %>% rbind(c(-1,0,0),c(0,0,1)) %>% seg %>% close %>% pl(fill= 'red')
  
  
  
  l(3,2,1) %>% pl
  l(3,2,0) %>% pl
  p(2,-3,1) %>% pl
  
  init()
  abline(h=0)
  abline(v=0)
  fta(ge, l(1,2,1), p(0,1,1)) %>% pl(col ='red')
  fta(ge, l(1,2,1), l(0,1,1)) %>% pl(col ='red')
  gell(shape = rmat(.8)) %>% pl -> ge0
  
  
  
  l(1,2,-1) %>% pl
  l(1,2,0) %>% pl
  
  fta(ge,l(1,2,-1), l(0,1,0)) %>% pl
  fta(ge0,l(1,2,-1), l(0,1,0)) %>% pl
  
  
  p(1,2,0) %>% pl
  ge %>% axes(p(1,2,0),'l',tan=T) %>% pl
  ge %>% axes(p(1,2,0),'l',tan=F) %>% pl
  ge %>% axes(p(3,1,1),'l',tan=T) %>% pl
  ge %>% axes(p(3,1,1),'l',tan=F) %>% pl
  
  ge %>% fta(l(0,1,0),p(0,-1,0)) %>% pl
  
  ge %>% fta(l(2,-1,1),p(1,1,0)) %>% pl(col = 'red')
  l(2,-1,1) %>% pl
  ge %>% axes(p(1,2,0)) %>% pl(col = 'red')
}


if(F){
  init()
  p1 <- p(2,3,1)  
  l1 <- l(1,2,1)  
  l2 <- l(2,-2,1)
  p1 %>% pl
  l1 %>% pl
  l2 %>% pl
  fta(p1,l1) %>% pl
  fta(p1,l1,l2) %>% pl %>% text(c('a','b'))
  fta(l1,p1,l2) %>% pl %>% text(c('a','b'),pos =3)
  fta(l1,p1,p(0,-1,0)) %>% pl %>% text(c('a','b'),pos =3)
  fta(l1,p1,p(0,1,0)) %>% pl %>% text(c('a','b'),pos =3)
  fta(l1,p1,p(0,-1,0)) 
  fta(l1,p1,p(0,1,0))
  fta(l1,l2,l(1,0,0)) %>% pl(col='red')
  gell() %>% pl

}

COND <- FALSE
if(COND) {
  
  # draw an ellipse
  
  init()

  gell(shape = cbind(c(2,.8),c(.8,1)), center=c(-2,1.5)) %>% pl ->ge
  
  # segment to x axis parallel to y axis
  
  ge %>% fta(l(0,1,0), p(0,-1,0)) %>% pl
  
  # segment to x axis with tangent parallel to x axis
  
  ge %>% fta(l(0,1,0), l(0,1,0)) %>% pl
  
  # segment to x axis with tangent parallel to y axis
  
  ge %>% fta(l(0,1,0), l(1,0,0)) %>% pl
  
  gell_box_(ge) %>% pl
  
  lightgray <- '#88888877'
    
  axes(ge,p(1,0,0),'l', tan = T) %>% pl
  axes(ge,p(1,0,0),'l', tan = F) %>% pl
  
  {
    # Drawing the canonical confidence ellipse
    # 
    init()
    abline(h=0,col=lightgray)
    abline(v=0,col=lightgray)
    gell(shape = diag(c(1.2,1))%*%rmat(-.8)%*%diag(c(1.2,1)), center = c(-.8,2.5)) %>% pl -> ge
    fta(ge, l(0,1,0), p(0,-1,0)) %>% pl %>% sel(2) %>% 
      text(lab('$\\hat{\\beta}_1$'),pos = 1, cex = 2)
    fta(ge, l(0,1,0), l(0,1,0)) %>% pl %>% sel(2) %>% 
      text(lab('$\\hat{\\gamma}_1$'),pos = 1, cex = 2)
    fta(ge, l(0,1,0), l(1,0,0)) %>% pl %>% sel(2) %>% 
      text(lab('$\\hat{\\psi}_1$'),pos = 1, cex = 2)
    l(0,1,-1) %>% pl
    fta(ge, l(0,1,-1), l(0,1,0) ) %>% pl(col = 'red') %>% 
      sel(2) %>% p %>%  fta(l(0,1,0)) %>% pl(col= 'red') %>% 
      sel(2) %>% p %>% text(lab('$\\hat{\\delta}_1$'), cex = 2, pos = 1)
    ge %>% gell_box_ %>% pl(col  = lightgray)
    fta(ge, l(1,0,0), p(1,0,0)) %>% pl(col = 'blue') -> vertax
    vertax %>% sel(2) %>% text(lab('$\\beta_2$'), cex = 2, pos=4, col = 'blue')
    fta(ge, l(1,0,0), l(1,0,0)) %>% pl(col = 'blue') -> vertax2
    vertax2 %>% sel(2) %>% text(lab('$\\gamma_2$'), cex = 2, pos=4, col = 'blue')
    }
  
  {
    # Drawing the acceptance ellipse and Simpson's Paradox
    # 
    inity()
    cex <-1.8
    abline(h=0,col=lightgray)
    abline(v=0,col=lightgray)
    gell(shape = diag(c(1.2,1))%*%rmat(-.8)%*%diag(c(1.2,1))) %>% pl -> ge
    
    gell_box_(ge) %>% pl(col=lightgray)
    
    fta(ge,l(0,0,-1), l(0,-1,0)) %>% rbind(c(0,1,0)) %>% seg %>%  close %>% pl(fill='#AAAAAA33')
    fta(ge,l(0,0,1), l(0,1,0)) %>% rbind(c(0,-1,0)) %>% seg %>%  close %>% pl(fill='#AAAAAA33')
    
    fta(ge,l(0,0,-1), l(-1,0,0)) %>% rbind(c(1,0,0)) %>% seg %>%  close %>% pl(fill='#AAAAAA33')
    fta(ge,l(0,0,1), l(1,0,0)) %>% rbind(c(-1,0,0)) %>% seg %>%  close %>% pl(fill='#AAAAAA33')
    
    text(p(-4,1,1),lab('$\\beta_2 > 0$'), pos = 3, cex = cex)
    text(p(-4,1,1),lab('$\\gamma_2 < 0$'), pos = 1, cex = cex)
    
    
    text(p(-1,3,1),lab('$\\beta_1 < 0$'), pos = 3,  cex = cex)
    text(p(-1,3,1),lab('$\\gamma_1 > 0$'), pos = 1, cex = cex)
    
    
    text(p(3,-1,1),lab('$\\beta_2 < 0$'), pos = 3, cex = cex)
    text(p(3,-1,1),lab('$\\gamma_2 > 0$'), pos = 1, cex = cex)
    
    text(p(1,-3,1),lab('$\\beta_1 > 0$'), pos = 3, cex = cex)
    text(p(1,-3,1),lab('$\\gamma_1 < 0$'), pos = 1, cex = cex)
  }
  
  
  {
    # Drawing the acceptance ellipse and Suppression Cone
    # 
    init()
    lightgray <- '#99999999'
      cex <-1.8
      abline(h=0,col=lightgray)
      abline(v=0,col=lightgray)
      gell(shape = diag(c(1.2,1))%*%rmat(-.4)%*%diag(c(1.2,1))) %>% pl -> ge
      gell(shape = diag(c(1.2,1))%*%diag(2)%*%diag(c(1.2,1))) %>% pl -> gc
      
      gell_box_(ge) %>% seg %>% close %>% pl(col=lightgray)
      
      
      lab <- latex2exp::TeX
      
      fta(ge,l(0,0,-1), l(0,-1,0)) %>% rbind(c(0,1,0)) %>% seg %>%  close %>% pl(fill='#AAAAAA33')
      fta(ge,l(0,0,1), l(0,1,0)) %>% rbind(c(0,-1,0)) %>% seg %>%  close %>% pl(fill='#AAAAAA33')
      
      fta(ge,l(0,0,-1), l(-1,0,0)) %>% rbind(c(1,0,0)) %>% seg %>%  close %>% pl(fill='#AAAAAA33')
      fta(ge,l(0,0,1), l(1,0,0)) %>% rbind(c(-1,0,0)) %>% seg %>%  close %>% pl(fill='#AAAAAA33')
      
      text(p(-3,.5),lab('$\\beta_2 > 0$'), pos = 3, cex = cex) %>% 
        text(lab('$\\gamma_2 < 0$'), pos = 1, cex = cex)
      
      
      text(p(-1/2,2,1),lab('$\\beta_1 < 0$'), pos = 3,  cex = cex) %>% 
        text(lab('$\\gamma_1 > 0$'), pos = 1, cex = cex)
      
      
      text(p(2.5,-1/2,1),lab('$\\beta_2 < 0$'), pos = 3, cex = cex) %>% 
        text(lab('$\\gamma_2 > 0$'), pos = 1, cex = cex)
      
      text(p(1/2,-2,1),lab('$\\beta_1 > 0$'), pos = 3, cex = cex) %>% 
        text(lab('$\\gamma_1 < 0$'), pos = 1, cex = cex)
      
      ints <- inter_gells(ge,gc)
      # fta(p(0,0), l(0,0,1), p(ints[1,])) %>% pl
      fta(p(0,0), l(0,0,1), p(ints[2,])) %>% 
        rbind(fta(p(0,0), l(0,0,1), p(ints[3,])) %>% flip) %>% seg %>% 
        pl(fill = '#88008811')
      fta(p(0,0), l(0,0,1), p(ints[4,])) %>% 
        rbind(fta(p(0,0), l(0,0,1), p(ints[1,])) %>% flip) %>% seg %>% 
        pl(fill = '#88008811')
      text(p(2,-2),'cone of SSR\nsuppression', cex = 2)
      text(p(-2,2),'cone of SSR\nsuppression', cex = 2)
  }
  
  
  {
    # Change score simulation
    
    mkd <- function() { 
      data.frame(devpre = rnorm(1000)) %>% 
        within({
          devpost <- .6*devpre + sqrt(1-.6^2) * rnorm(1000)
          X1 <- rep(c(0,1), each = 1000/2)
          X2 <- 5 + 2*X1 + devpre
          Y <- 5 + 2*X1 + devpost
        })
    }
    dd <- mkd()
    dd <- sortdf(dd, ~ X2)
    library(latticeExtra)
    xyplot(Y ~ X2, dd, groups = X1)
    fit <- lm(Y ~ X1 + X2, dd)
    dd $ pred <- predict(fit)
    xyplot(Y ~ X2, dd, groups = X1) +
      xyplot(pred ~ X2, dd, groups = X1, type = 'l')
    inite(xlim = c(-1,2))
    abline(h=0)
    abline(v=0)
    cell(fit) %>% lines
    ge <- gell(cell(fit))
    init()
    
    ge %>% pl(col='red') 
    l(0,-1,1) %>% pl
    fta(ge, l(0,-1,1), l(0,1,0)) %>% pl %>% sel(2) -> del1
    fta(p(del1),l(0,1,0)) %>% pl
    fta(p(ge$center), l(0,1,0)) %>% pl
    l(1,1,-1) %>% rot(pi/3) %>% pl %>% fta(p(-2,2), .) %>% pl
    rbind(c(-1,-1,1),c(2,-1,0), c(5,1,0), c(1,1,1)) %>% 
      seg %>% close %>% pl(fill= '#00777722')
    rbind(c(-1,-1,1),c(2,-1,0), c(5,1,0), c(3,1,1)) %>% 
      seg %>% close %>% pl(fill= '#77007722')
    
  }
  
  {
    ## Combining estimators and line of osculation ----
    ## 
    ## 
    init(xlim = c(-2,6))
    library(spida2)
    gell(shape = rmat(.8)) %>% pl -> g1
    gell(shape = rmat(-.7), center = c(4,-1)) %>% pl -> g2
    
    
    comb(g1,g2) %>% pl
    fta(g1,g2) %>%  pl
    
    # How evidence overwhelms a prior
    
    init(xlim=c(-2,5))
    gell(radius=3) %>% pl(col = 'red')->gprior
    lab <- latex2exp::TeX
    p(2,-2,1) %>% text(lab('prior'), col = 'red')
    
    freq_ell <- lapply(c(2,5,10,25,50,100),
                       function(n) {
                         gell(shape=rmat(-.9), center = c(3,4), radius = 1/sqrt(n))
                       }
    )
    freq_ell %>% pl
    post_ell <- lapply(freq_ell, function(ge) comb(ge,gprior))
    post_ell %>% pl(col='green')
    
    
    lg <- '#88888833'  
    
    
    for(lam in seq(0,1,.01)) {comb(g1,g2,lam) %>% {c(.$center)} %>%  p %>% pl}
    for(lam in seq(0,1,.01)) {comb(g1,g2,lam) %>%  pl(col = '#77777722')}
    comb(g1,g2)    %>% pl (fill = lg)
    fta(g1,g2) %>%  pl(col = hcl(0,50,100,.5))
  }
  
  {
    # playing with colours
    pal
    
    se <- seq(0,150,10)/150
    hcl(seq(0,360,10), 100,100) %>% pal
    hcl(240, se,50) %>% pal
    
    hsv(.5, .5,.5) %>% pal
    hsv(se, .5, 1) %>% pal
    hsv(.5, se, 1) %>% pal
    hsv(.5, 1,se) %>% pal
    
    
    
    lightgray %>% col2rgb
  }
  
  {
    # Paradoxical shrinkage
    init()
    gell(A=rotm(pi/20)%*%cbind(c(-1,2),3*c(.1,.1))) %>% pl -> gprior
    gell(shape=rmat(-.98), center = c(2.5,-.2)) %>% pl -> gdata
    fta(gprior, gdata) %>% pl
    wei <- 0
    {
      wei <- wei + .1
     comb(gprior, gdata, wei) %>% pl(col = 'red') %>% cen  %>% pl(pch = 16, col = 'red')
    }
    for(rad in seq(.1,5,.1)){
      gell(A=rotm(-pi/8)%*%cbind(c(-1,2),3*c(.1,.1)), radius = rad) %>% pl(col=lightgray) -> gprior
      gell(shape=rmat(-.98), center = c(2.5,-.2), radius = rad) %>% pl(col=lightgray) -> gdata
      
    }
  
    comb(gprior, gdata, .3) %>% pl %>% cen %>% pl(pch =15, col ='red')
    
    
  }
  
  {
    # intersecting a function with a value
    # 
    init()
    fun <- function(x,y) x^3 - 2*y^3 + x*y
    fun <- function(x,y) sin(x)^3 - sin(y)^3
    inter(fun,0) %>% 
      Map(pl,.,fill = c('red','blue','green'))
  }
  
}





