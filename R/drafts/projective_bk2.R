#'
#' THis is a part of a package for plotting objects
#' represented in projective 2 dimensional space
#' Using objects 'ppoint' and 'pline'
#' 
#' General approach
#' ---
#' title: Projective approach to Generalized Ellipsoids
#' ---
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
#'   orthogonal $\Gamma$ and singular values which may be
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
#' - gell and
#'   - point
#'     - scale to gell
#'     - conjugate axis
#'     - conjugate parallelogram
#'     - subtending rectangle
#' - gell and
#'   - line
#'     - tangent point
#' - gell and gell
#'   - locus of osculation
#'        
#'     
#'     
#' 
#'

library(p3d)
library(spida2)
library(latex2exp)

## init ###################################
init <- function(from = -fac*c(1,1), to = fac*c(1,1), fac = 3,
                 xlab = '', ylab = '', asp = 1) {
  plot(rbind(from,to), type ='n', xlab = xlab, ylab = ylab, asp = asp)
}
COLS <- c('red','cyan','purple','yellow','teal' ,'blue','#4444FF','green','goldenrod',
               '#23423411',
               '#014488')
               
rss <- function(x) {
  # root sum of squares: vector length
  M <- max(abs(x))
  if(M == 0) return(0)
  M * sqrt(sum((x/M)^2))
}


n <- function(x) {
  if(ncol(x) != 3) stop('can only normalize 3-vectors')
  div <- ifelse(x[,3] == 0, pracma::hypot(x[,1],x[,2]), x[,3])
  x/div
}

rbind(c(1,1,0),c(1,0,1),c(0,0,0)) %>% n

pl <- function(x,...) {
  # plot method
  UseMethod('pl')
}

r <- function(...) {
  rbind(c(...))
}

##+ line ####


l <- function(x,...) {
  # constructor for 'line'
  UseMethod('l')
} 


l.line <- function(x,...) {
  x
}

l.default <- function(x,...) {
  # 
  # single line: l(1,1,0) or l(c(1,1,0))
  # multiple: l(c(1,1,0),c(1,1,1))
  # or: l(rbind(c(1,1,0),c(1,1,1)))
  # 
  if(length(x) == 1) x <- rbind(c(x, ...))
  else x <- rbind(x,...)
  if(ncol(x) == 2) x <- cbind(x,0)
  x <- n(x)
  class(x) <- 'line'
  x
}

setOldClass('line')


pl.line <- function(x,...) {
  x <- l(x)
  for(i in seq_len(nrow(x))) {
    ll <- x[i,]
    if(ll[2] != 0) abline(a = -ll[3]/ll[2], b = -ll[1]/ll[2], ...)
    else abline(v = -ll[3]/ll[1], ...)
  }
  invisible(x)
}

init()
l(1,1,0) %>% pl

##+ point constructors ####



p <- function(x,...) {
  # constructor
  UseMethod('p')
} 

p.default <- function(x,...) {
  # 
  # single line: l(1,1,0) or l(c(1,1,0))
  # multiple: l(c(1,1,0),c(1,1,1))
  # or: l(rbind(c(1,1,0),c(1,1,1)))
  # 
  
  if(length(x) == 1) x <- rbind(c(x, ...))
  else x <- rbind(x,...)
  if(ncol(x) == 2) x <- cbind(x,1)
  x <- n(x)
  class(x) <- 'point'
  x
}

setOldClass('point')

pl.point <- function(x,...) {
  x <- p(x)
  #disp(x)
  x <- x/x[,3]
  points(x[,1:2,drop=FALSE], ...)
  # for(i in seq_len(nrow(x))) {
  #   ll <- x[i,,drop = FALSE]
  #   if(ll[3] == 1) points(ll[,1:2,drop=FALSE],...)
  # }
  invisible(x)
}

##' join #### 

# replaces 'dual' for some purposes

setGeneric('join', function(x, y, ...) standardGeneric('join'))

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

## segments ####

seg <- function(x, ...) {
  #  constructor for segments
  UseMethod('seg')
}

setOldClass('seg')

seg.line <- function(x, ...) {
  if(nrow(x) <3) stop('need 3+ lines')
  x <- rbind(x, x[1,])
  ret <- matrix(0, 0, 3)
  for(i in 2:nrow(x)) ret <- rbind(ret, join(x[i-1,],x[i,]))
  class(ret) <- 'seg'
  ret
}

seg.point <- function(x, ..., close = FALSE) {
  if(nrow(x) <2) stop('need 2+ points')
  if(close) x <- rbind(x, x[1,])
  ret <- x
  class(ret) <- 'seg'
  ret
}

seg.seg <- function(x, ...) {
  x
}

seg.default <- function(x, ...) {
  x <- p(x)
  class(x) <- 'seg'
  x
}


in_seg <- function(X, seg2) {
  # X is a set of points, usually on the line spanning seg2
  # seg2 is a 2-point segment
  # which points in X have projections onto line
  # spanning seg2 that are within seq2 
  # 
  # Can be dually defined to identify lines in cones
  X <- p(X)
  seg2 <- seg(seg2)
  line <- l(seg2)
  u_perp <- l(c(-line[1,2],line[1,1], 0))  # perpendicular
  Xpos <- X %*% t(u_perp)
  Xpos <- Xpos/X[,3]
  segpos <- seg2 %*% t(u_perp)
  segpos <- segpos/seg2[,3]
  low <- min(segpos)
  high <- max(segpos)
  (Xpos >= low) & (Xpos <= high)
}

order_seg <- function(X, seg2 = X[c(1,2),]) {
  # X is a set of points, usually on the line spanning seg2
  # seg2 is a 2-point segment
  # which points in X have projections onto line
  # spanning seg2 that are within seq2 
  # 
  # Can be dually defined to identify lines in cones
  X <- p(X)
  seg2 <- seg(seg2)
  line <- join(p(seg2[1,]),p(seg2[2,]))
  u_perp <- l(c(-line[1,2],line[1,1], 0))  # perpendicular
  Xpos <- X %*% t(u_perp)
  Xtheta <- atan2(Xpos, X[,3])
  segpos <- seg2 %*% t(u_perp)
  segtheta <- atan2(segpos, seg2[,3])
  # low <- min(segtheta)
  # high <- max(segtheta)
  (Xtheta -segtheta[1]) / (segtheta[2] - segtheta[1] )
}
init()
seg2 <- seg(p(rbind(p(1,1,1),p(3,2,1))))
rbind(
  c(1,1,1),
  c(3,2,1),
  c(2,1,0),
  c(21,11,1),
  c(2001,1001,1),
  c(-2,-1,0)
) %>% p %>% seg %>% order_seg(seg2)

finite.seg <- function(x, radius = 10, n =30) {
  # 
  # points at infinity are replaced with finite points on a circle of radius
  # assumes that first and last are finite points
  # assumes that circle includes all finite points
  # 
  # If a point is at infinity while the previous one was not,
  # replace it with the corresponding point on the circle.
  # 
  # If the next point is finite, add a point on the circle
  # from the point at infinity
  # 
  # Summary:
  # 
  # For each point:
  # - if finite add
  # - if infinite: 
  #   - check previous
  #     - if finite add point on circle (need to track)
  #     - if infinite, ignore
  #   - check next:
  #     - if finite: add point on circle (track)
  #     - if infinite, ignore
  # - go to next point
  # - Reloop:
  # - fill in points between points on circle
  # 
  # Checks:
  # - first and last should be finite
  x <- seg(x)
  el <- gell.default(radius = radius)
  if(x[1,3] == 0 || x[nrow(x),3] == 0) stop('First and last points should be finite')
  if(all(x[,3] != 0)) return(x)
  ret <- NULL
  for(i in 1:nrow(x)) {
    if(x[i,3]!=0) {
      ret <- rbind(ret, x[i,,drop=FALSE])
    } else {  # point at infinity
      # check previous
      if(x[i-1,3] != 0) {
        ret <- rbind(ret, 
                     gell_intersects(el, join(p(x[i-1,]),p(x[i,])))[2,])
      }
      # check next  
      if(x[i+1,3] != 0) {
          ret <- rbind(ret, 
                       gell_intersects(el, join(p(x[i+1,]),p(x[i,])))[2,])
      }
    }
  }
  # add circumferential points
  seg(ret)
}

rbind(
  c(1,1,1),
  c(1,2,1),
  c(1,1,0),
  c(3,1,1)
) %>% seg %>% finite.seg



pl.seg <- function(x, ..., type = 'l') {
  # NOTE: Uses gell_intersects below
  # finite segments are easy
  # challenge is drawing segments that
  # involve points at infinity
  # We do it here by defining a circle 
  # larger than the viewport and 
  # drawing segments that are within the circle
  # and, for points at infinity, drawing the
  # segment to the point on the circle on
  # the line that points to the point at infinity
  # 
  if(nrow(x) < 2) stop('need 2+ rows')
  
  usr <- par('usr')
  bigell <- gell.default(radius=sum(abs(usr)))
  
  for(i in 1:(nrow(x) - 1)) {
    p1 <- p(x[i,])
    p2 <- p(x[i+1,])
    seg <- p(x[i:(i+1),,drop=FALSE])
    if(any(seg[,3] == 0)) {
      if(all(seg[,3] == 0)) next   # line at infinity
      ord <- (seg[,3] == 0) + 1  # to order finite point first
      seg[] <- seg[ord,] 
      lin <- join(p(seg[1,]),p(seg[2,]))
      circpts <- gell_intersects(bigell, lin)
      segs <- rbind(seg, circpts)
      disp(segs)
      os <- order_seg(segs)
      keep <- os >= 0
      keep[2] <- FALSE
      segs <- segs[keep,,drop =FALSE]
      lines(segs[,-3], type = type , ...)

      ## %>% HERE %>% #####
      
    } else {
      lines(seg[,-3], type = type, ...)
    }
  }
}
init()
seg(p(rbind(p(1,1,1),p(3,-2,1),p(-3,3,1)))) %>% pl(type='b')
seg(p(rbind(p(1,1,1),p(1,1,0),p(-1,-1,1)))) %>% order_seg
seg(p(rbind(p(1,1,1),p(1,0,0),p(-1,-1,1))), close = T) %>% mult(rot(-2*pi/100)) %>% pl
seg(p(rbind(p(1,1,1),p(1,0,0),p(-1,-1,1))), close = T) %>% mult(rot(1*pi/100),.) %>% pl
gell.default() %>% plus(p(c(2,2))) %>% pl
gell.default() %>% plus(p(c(2,2))) %>% mult(rot(-2*pi/100)) %>%pl
p(rbind(p(1,1,1),p(3,4,1))) %>% join

rbind(c(1,1,1),c(1,0,0),c(2,-1,1),c(1,1,1)) %>% seg %>% pl.seg


# seg(p(rbind(p(1,1,1),p(3,4,1),p(-3,3,1)))) %>% p %>% join %>% pl(pch=16)

if(FALSE) {
  
  # test
  seg(p(p(1,1,0),p(1,1,1), p(2,1,1))) %>% l
  
  
  
  init()
  p(1,1) %>% pl
  p(1,2) %>% pl
  join(p(1,1),p(1,2)) %>% pl(col='red')
  join(p(c(1,1),c(1,-1)), p(c(2,2))) %>% pl
  join(
    p(c(1,1,0)),  # point at infinity)
    p(c(1,2))
  ) %>% pl
  join(p(c(1,1,1)), p(c(1,1,1))) %>% join(l(c(1,-1,0))) %>% join(p(c(1,1,1))) %>% pl
}

# p.line <- function(l1, l2, ...) {
#   # redundant should be same as join
#   # intersection(s) of l1 and l2 
#   if(missing(l2)) {
#     if(nrow(l1) == 2) p(l(l1[1,]),l(l1[2,]))
#     else stop('need 2 lines')
#   } else {
#     dual(l1,l2)
#   }
# }

# l.point <- function(p1, p2,...) {
#   # constructor for line from two points
#   # redundant: should be same as join
#   if(missing(p2)) {
#     if(nrow(p1) == 2) l(p(p1[1,]),p(p1[2,]))
#     else l(l(p1,p(c(0,0,1))))
#   } else {
#     dual(p1,p2)
#   }
# }




# test
plot(c(-5,5),c(-5,5), type = 'n', asp = 1)

p(c(0,0)) %>% 
  join(p(c(1,1))) %>% pl -> line1
p(c(1,1)) %>% 
  join(p(c(-1,1))) %>% pl -> line2
line1 %>% join(line2) %>% pl(pch=16, col = 'red')

##' DO LATER ####

pl.ell <-function(x, ...) {
  lines(x,...)
  invisible(x)
}

seg <- function(x, ...) {
  UseMethod('seg')
}

seg.line <- function(ll, ...) {
  ll <- l(ll)
  if(nrow(ll) < 2) stop('need at least two lines')
  ind <- seq_len(nrow(ll))
  # disp(ind)
  pts <- dual(ll,ll[c(ind[-1],ind[1]),, drop = FALSE]) %>% p
  pts <- rbind(pts,pts[1,]) %>% p
  class(pts) <- 'seg'
  pts
}


pl.seg <- function(x, ...) {
  # 
  # need to deal with points at infinity better
  # 
  for(i in seq_len(nrow(x))) {
    if(x[i,3] == 0) x[i,1:2] <- NA
  }
  lines(x[,-3,drop = FALSE], ...)
}

##' END DO LATER ####

## plus, mult ####


Ops.point <- function(e1, e2) {
  if(.Generic == '+') {
    return(plus(e1,e2))
  }
  e1 <- unclass(e1)
  e1 + e2
} 

p(0,1,1) %>% plus( p(1,3,1))

 
setOldClass('matrix')
setGeneric('plus', function(x,y,...) standardGeneric('plus') )
setGeneric('mult', function(x, y, ...) standardGeneric('mult') )

setMethod('mult', signature('line','numeric'),
          function(x, y, ...) {
            x <- l(x)
            l(cbind(x[,1:2,drop=FALSE]/y,x[,3]))
          } 
)
setMethod('mult', signature('numeric','line'),
          function(x, y, ...) {
            mult(y, x, ...)
          } 
)

setMethod('mult', signature('line','matrix'),
          function(x, y, ...) {
            x <- l(x)
            l(cbind(x[,1:2,drop=FALSE]%*%t(solve(y)),x[,3]))
          } 
)

setMethod('mult', signature('matrix','line'),
          function(x, y, ...) mult(y,t(x))
)




setOldClass('seg')

setMethod('mult', signature('point','matrix'),
          function(x, y, ...) {
            x <- p(x)
            p(cbind(x[,1:2,drop=FALSE] %*% y, x[,3]))
          } 
)
setMethod('mult', signature('matrix','point'),
          function(x, y, ...) mult(y,t(x),...)
)
setMethod('mult', signature('seg','matrix'),
          function(x, y, ...) {
            x <- seg(x)
            seg(cbind(x[,1:2,drop=FALSE] %*% y, x[,3]))
          } 
)
setMethod('mult', signature('matrix','seg'),
          function(x, y, ...) mult(y,t(x),...)
)



rot <- function(theta) {
  cbind(c(cos(theta), -sin(theta)), c(sin(theta),cos(theta)))
}

if(FALSE) {
  rot(pi)
  init()
  rbind(
    p(2,1,1),
    p(1,2,1),
    p(-1,0,1)) %>% p %>% pl(pch=16) ->p1
  p1 %>% mult(rot(pi/2)) %>% pl(pch=16, col ='red') -> pr1
  join(p1,pr1) %>% pl
  
  mult(p(1,2,1),-2*diag(2))
  
  line2 %>% mult(2) %>% pl(col = 'red')
  mult(2, line2) %>% pl(col='blue')
  
  l(c(1,1,1))         
  l(c(1,1,1)) %>% mult(2) %>% pl
  
  l(c(1,1,0)) %>% pl(lwd=2,col = 'green')
  l(c(1,1,0)) %>% mult(rot(pi/2)) %>% pl(lwd=2,lty=2,col = 'green')
  for(i in 0:1000) {
    l(c(1,1,1)) %>% mult(rot(i*2*pi/1000)) %>% mult(cbind(c(2,1),c(-1,1))) %>% pl
  }
}

setMethod('mult', signature('point','numeric'),
          function(x, y, ...) {
            x <- p(x)
            p(cbind(x[,1:2, drop=FALSE]*y,x[,3]))
          } 
)
setMethod('mult', signature('numeric','point'),
          function(x, y, ...) {
            mult(y, x, ...)
          } 
)

p(c(-2,-5)) %>% pl %>% join(p(c(0,0))) %>% pl


setMethod(
  'plus', 
  signature('line','point'),
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

setMethod(
  'plus', 
  signature('point','line'),
  function(x, y, ...) {
    mult(y, x, ...)
  } 
)

setMethod(
  'plus', signature('point','point'),
  function(x, y, ...) {
    obj <- p(x)
    pt <- p(y)
    nr <- nrow(obj)
    ind <- seq_len(nrow(pt))
    ind <- rep(ind, length.out = nr)
    pt <- pt[ind,,drop=FALSE]
    disp(obj)
    disp(pt)
    ret <- obj[1:2,] * pt[,3] + pt[1:2,] * obj[,3]
    ret <- cbind(ret, pt[,3]*obj[,3])
    p(ret)
  } 
)

.plus <- function(x,y,...) {
  obj <- p(x)
  pt <- p(y)
  nr <- nrow(obj)
  ind <- seq_len(nrow(pt))
  ind <- rep(ind, length.out = nr)
  pt <- pt[ind,,drop=FALSE]
  disp(obj)
  disp(pt)
  ret <- obj[,1:2,drop=FALSE] * pt[,3] + pt[,1:2,drop=FALSE] * obj[,3]
  ret <- cbind(ret, pt[,3]*obj[,3])
  p(ret)
  
}

p(rbind(c(1,1,1), c(1,2,0), c(2,1,1), c(4,2,0))) %>% .plus( 
p(rbind(c(3,1,1), c(1,2,1), c(5,2,0), c(1,2,0)))) 




init()
c(1,1) %>% p %>% join(c(2,3) %>% p) %>% pl -> line1
c(0,0) %>% p %>% pl
line1 %>% plus(p(c(2,-1))) %>% pl(col='red')

line2 %>% mult(2) %>% pl(col = 'red')
mult(2, line2) %>% pl(col='blue')

p(c(0,0),c(1,1),c(-2,1)) %>% pl(col = 'green', pch = 16)
p(c(0,0),c(1,1),c(-2,1)) %>% plus(p(c(-2,-2))) %>% pl(col = 'red', pch = 16)



l(c(1,1,1))         
l(c(1,1,1)) %>% mult(2) %>% pl  



p(c(-2,-5)) %>% pl %>% join(p(c(0,0))) %>% pl

perp <- function(line, point = p(c(0,0))) {
  # line perpendicular through a given point
  line <- l(line)
  point <- p(point)
  ret <- line[,c(2,1,3), drop = FALSE]
  ret[,1] <- -ret[,1]
  ret[,3] <- 0
  plus(l(ret), point)
} 

if(FALSE){
  init()
  join(p(-1,-1),p(-2,2)) %>% pl -> line1
  p(p(-1,-1),p(-2,2)) %>% pl
  perp(line1) %>% pl(col = 'red')
  p(0,0) %>% pl(pch=16)
  perp(line1,p(2,2)) %>% pl(col = 'red')
}


rows <- function(m) {
  ret <- lapply(seq_len(nrow(m)), function(i) {
    ret <- m[i,,drop=FALSE]
    class(ret) <- class(m)
    ret
  })
  ret
}


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

dual.ell <- function(x, ...) {
  shape <- attr(x,'params')$shape
  radius <- attr(x,'params')$radius
  ell(shape = solve(shape), radius = 1/radius)
}

circle <- function(n = 100) {
  theta <- seq(0, 2*pi, length.out = n)
  cbind(cos(theta), sin(theta))
}

init()
circle() %>% l %>% mult(2) %>% pl
circle() %>% p %>% mult(2) %>% pl(lwd = 3, type = 'l', col = 'red')

text.point <- function(x, labels ='', ...) {
  x <- p(x)
  nr <- nrow(x)
  ret <- x[,1:2, drop = FALSE]
  keep <- x[,3] != 0
  ret <- ret[keep,, drop = FALSE]
  if(length(labels) > 1) labels <- rep(labels, length.out = nr)[keep]
  text(ret, labels = labels, ...)
  invisible(ret)
}

text(p(c(-3,0)), cex = 2,labels = TeX('here: $\\int_0^\\infty \\alpha_2 dx =\\beta \\times \\frac{x^2}{1+x^3}$'))

#'
##' ## Generalized ellipses ####
#'

gell <- function(x,...) {
  UseMethod('gell')
}
# gell.ell <- function(ell,...) {
#   params <- attr(m,'params')
#   ei <- eigen(params$shape)
#   ret <- list(gamma = ei$vectors, d = sqrt(ei$values), radius = params$radius, center = params$center)
#   class(ret) <- 'gell'
#   ret
# }

gell.default <- function(center = c(0,0), shape = diag(2), radius = 1, 
             gamma = eigen(shape)$vectors, d = sqrt(eigen(shape)$values)) {
  ret <- list(gamma = gamma, d = d, radius = radius, center = center)
  class(ret) <- 'gell'
  ret
}

gell.matrix <- function(X) {
  # data ellipse
  center <- apply(as.matrix(X), 2, mean)
  X <- t(t(X) - center)/sqrt(nrow(X))
  sv <- svd(X, nu = 0)
  ret <- list(gamma = sv$v, d = sv$d, radius = 1, center = center )
  class(ret) <- 'gell'
  ret
}  

pl.gell <- function(gell, ...) {
  if(max(gell$d) < Inf) {
    U <- circle()
    spida2::disp(gell$radius)
    spida2::disp(gell$d)
    spida2::disp(gell$center)
    pts <- t(gell$gamma %*% (t(U) * gell$radius * gell$d ) + gell$center)
    lines(pts, ...)
  }
  else {
    warning('pl.gell: case not yet implemented')
  }
  invisible(gell)
}

sell <- gell.default(shape = 0*diag(2) + 1, center = c(1,1))
sell %>% gell_intersects(seg(rbind(c(2,1,1),c(-2,2,1))))
sell %>% pl

setOldClass('gell')

setMethod(
  'plus', signature('gell','point'),
  function(x, y, ...) {
    ret <- x
    ret$center <- plus(p(ret$center), y)[1:2]
    ret
  } 
)
disp <- spida2::disp
setMethod(
  'mult', signature('gell','matrix'),
  function(x, y, ...) {
    ret <- x
    # print(ret)
    ret$center <- (x$center %>% p %>% mult(y))[,1:2] 
    #  disp(ret$center)
    sv <- svd(t(y) %*% ret$gamma %*% diag(ret$d), nv =0)
    ret$gamma <- sv$u
    ret$d <- sv$d
    ret
  } 
)
setMethod(
  'mult', signature('matrix','gell'),
  function(x, y, ...) mult(y, t(x), ...)
)

center <- function(gell,...) {
  # generic in spida2
  p(gell$center)
}
pcenter <- center
if(F){
  init()
  gell(X) %>% pl
  gell(X) %>% plus(p(1,2)) %>% mult(rot(0))%>% pl
}
X <- cbind(rnorm(200),rnorm(200)) %*% cbind(c(1,1), c(.2,-.3))
init()
X %>% p %>% pl  
gell(X) %>% pl
methods(pl)

gell(X) %>% center %>% pl(pch = 16)

zero <- function(x,...) {
  UseMethod('zero')
}
zero.gell <- function(g, ...){
  ret <- g
  ret$center <- c(0,0)
  ret
}


## axes ####

setGeneric('axes', function(x,y,...) standardGeneric('axes'))

setMethod(
  "axes",
  signature(x = 'gell', y = 'point'),
  function(x, y, ...)  {
    P <- cbind(c(0,1),c(-1,0))
    pt <- p(y)
    gell <- x
    pt <- plus(pt, mult(pcenter(gell),-1)) # centered ellipse
    v <- diag(1/gell$d) %*% t(gell$gamma) %*% cbind(pt[1:2])
    u <- v / sqrt(sum(v^2))
    pts <- gell$gamma %*% diag(gell$d) %*% cbind(u, P%*%u)
    pts_on_ell <- p(t(pts)) %>% plus(center(gell))
    tan_dirs <- p(t(pts))[2:1,]
    tan_dirs[,3] <- 0
    list(points = pts_on_ell, tan_dirs = p(tan_dirs))
  })

sel <- function(x, nam) x[[nam]]

setOldClass('gell')
setOldClass('point')
setOldClass('line')

X <- t(c(1,2)+t(X))
init()
gell(X) %>% pl
gell(X) %>% zero %>% pl

gell(X) %>% axes(p(1,1)) %>% sel('points') %>% pl(col = 'red')
gell(X) %>% axes(p(1,1)) %>% sel('tan_dirs') %>% join(center(gell(X))) %>% pl
gell(X) %>% center %>% pl(col = 'red')
gell(X) %>% axes(p(1,0,0)) %>% lapply(pl)
gell(X) %>% axes(p(1,2,1)) %>% lapply(pl)

# 
# axes yields the intersection of the ellipse with a segment from the
# center of the ellipse towards a point, at infty or finite
# - also the point on the ellipse that is conjugate
# 
# 

circle_intersects_1 <- function(y) {
  y <- l(y)
  e <- sum(y[1,1:2]^2 - 1)
  if(y[1,3] != 0 && e < 0) return(p(0,0,0))
  xs <- (-c(1,-1)*y[1,2]*sqrt(e) - y[1,1])/(e+1)
  ys <- -(y[1,1]*xs +1)/y[1,2]
  p(cbind(xs,ys,1))
}

circle_intersects <- function(y) {
  y <- l(y)
  # if line horizontal
  l2 <- rss(y[1,1:2])
  if(y[1,3] != 0 && l2 < 1) return(p(matrix(0,0,3)))
  theta1 <- asin(-y[1,3]/l2)
  pts1 <- p(cbind(c(1,-1)*cos(theta1),sin(theta1),1))
  pts <- pts1 %>% mult(rot(atan2(-y[1,1],y[1,2])))
  pts
}
 
gell_intersects <- function(e, obj) {
  if(class(obj) == 'line') invisible(gell_intersects_line(e, obj))
  else if(class(obj) == 'seg') invisible(gell_intersects_seg(e, obj))
}

gell_intersects_line <- function(e, line) {
  # seems to work !!!
  cen <- e$center
  ez <- e
  ez$center <- c(0,0)
  line <- line %>% plus(p(-cen))
  circle_line <- line %>% mult((ez$gamma)) %>% mult(diag(1/ez$d)) %>% mult(1/ez$radius)
  circle_pts <- circle_intersects(circle_line)
  pts <- circle_pts %>% mult(diag(ez$d)) %>% mult(t(ez$gamma)) %>% mult(ez$radius) %>% plus(p(cen))
  pts
}
gell_intersects_seg <- function(e, sg) {
  sg <- seg(sg)
  if(nrow(sg) != 2) stop('need segment with two rows')
  pts <- gell_intersects_line(e, join(p(sg[1,]),p(sg[2,])))
  sgs <- rbind(sg,pts)
  os <- order_seg(sgs)
  keep <- os >= 0 & os <= 1
  keep[1:2] <- FALSE
  p(sgs[keep,,drop=FALSE])
}

init()
gell.default(radius = 3) %>% pl
gell.default(radius = 3) %>% gell_intersects(l(1,1,1)) %>% pl
l(1,1,1) %>% pl
l(1,1,1) %>% join(l(1,0,0)) %>% pl
gell.default(radius = 3) %>% gell_intersects(seg(rbind(c(1,1,1),c(1,2,1)))) %>% pl
gell.default(radius = 3) %>% gell_intersects(seg(rbind(c(1,1,1),c(1,5,1)))) %>% pl
seg(rbind(c(1,1,1),c(1,1,0),c(-1,-1,1))) %>% pl




init()
el <- gell.default(shape = .1*diag(2) + .3, radius = 3, center = c(1,1))

el %>% pl
l(-3,3,1) %>% pl
el %>% gell_intersects(l(-3,3,1)) %>% pl
el %>% gell_intersects(l(-.1,.1,1)) %>% pl
l(-.1,.1,.01) %>% gell_intersects(el,.) %>% pl
seg(rbind( c(1,1,1),c(3,4,1))) %>% pl
seg(rbind( c(1,1,1),c(3,4,1))) %>% gell_intersects(el,.) %>% pl
seg(rbind(c(1,1,1), c(1,0,0))) %>% pl
seg(rbind(c(1,1,1), c(1,0,0))) %>% gell_intersects(el,.) %>% pl
seg(rbind(c(1,1,1), c(-1,0,0))) %>% gell_intersects(el,.) %>% pl

-p(0,0,1) %>% pl
init()
l(1,0,0) %>% pl
l(0,1,0) %>% pl
gell() %>% pl
l(-2,-2,1) %>% pl
l(-3,-2,1) %>% pl
circle_intersects(l(-3,-2,1)) %>% pl %>% print

l(1,1,0) %>% pl
l(1,1,0) %>% circle_intersects %>% pl
l(1,0,0) %>% circle_intersects %>% pl





l(-1,0,1) %>% pl
l(0,1,1) %>% pl
l(0,-1,1) %>% pl




gell_conj <- function(ge, dir) {
  # line thru center conjugate to a direction
  # i.e. intersecting ellipse with tangent parallel to a direction
  a <- axes(ge, dir)
}

gell_box <- function(ge, p1, p2) {
  a1 <- axes(ge, p1)
  if(missing(p2)) {
    a2 <- p(a1[c(2,1),])
  } else {
    a2 <- axes(ge, p2)
  } 
  l1 <- join(p(a1[1,]), p(a1[1,]) %>% plus)
  
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

plus.gell <- function(gell, pt)

plot(X)
gell(X) %>% pl 
gell(X) %>% axes(p(-c(0,1))) %>% pl(pch = 16, col = 'red')
gell(X) %>% axes(-c(0,1)) %>% pl(pch = 16, col = 'red')

# let's see whether points at infinity work as they should
# 
X <- t(c(2,3) + t(X))

X %>% plot
gell(X) %>% pl(col = 'red')
gell(X) %>% axes(c(1,0)) %>% pl(pch = 16, col = 'red')
gell(X) %>% pcenter %>% pl(pch = 16, col = 'blue')
gell(X) %>% axes(p(c(1,0,0))) %>% pl(pch = 16, col = 'red')


dual.gell <- function(gell, center = gell$center) {
  ret <- gell
  ret$d <- 1/ret$d
  ret$center <- center
  ret
}

X %>% plot(asp = 1)
gell(X) %>% pl(col = 'red')

gell(X) %>% dual %>% pl
(center = pcenter(gell(X))[[1:2]])

#############  HERE  #################################

center.ge

pl(gell(X))

svd(diag(c(3,4)))



## init ###################################
init <- function(from = -fac*c(1,1), to = fac*c(1,1), fac = 3,
                 xlab = '', ylab = '', asp = 1) {
  plot(rbind(from,to), type ='n', xlab = xlab, ylab = ylab, asp = asp)
}
COLS <- c('red','cyan','purple','yellow','teal' ,'blue','#4444FF','green','goldenrod',
               '#23423411',
               '#014488')


Sigma <- rbind(c(1.1, .9), c(.9, 1.0))
zell <- ell(c(1,1), shape = Sigma, radius = 1)

{
  
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


{
  
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


}
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

               

               
        






init();abline(v=-10:10,col='gray50',lty=2);abline(h=-10:10,col='gray50',lty=2)

zp <- p(rbind(c(0,1),c(1,1), c(-2,2)))
Map(pl,rows(zp),col=1:3, pch = 16, cex = as.list(1:3))
Map(pl,rows(zp %>% plus(p(c(-2,-2)))),col=1:3, pch = 16, cex = as.list(1:3))

zp %>% pl( pch = 16, col = 1:3)
zp 
zp %>% apply(1, print)

zp %>% plus(p(c(1,1,0)))
zp %>% pl

zl <- l(c(1,1,0),c(1,1,2), c(.5,2,1))
zl %>% pl
zl %>% plus(p(c(2,2))) %>% pl(col='red')

dual(zl[1,],zl[3,]) %>% p %>% pl(cex=5, col = 'orange', pch = 16)
dual(zp,zp[c(2,3,1),]) %>% l %>% pl




abline(v=Inf)


pl.default <- function(x,...) {
  
}

      if(ll[2] == 0) return(NULL)
      else abline(v=-ll[3],...)
    } else {
      if(ll[2] == 0) abline(h = -ll[3],...)
      else abline()
    }
      
    }
  }
}

#' - 
# library(p3d)
library(pracma)
library(spida2)
?ellplus
#' Drop 3rd coordinate of 3-vector(s)
#'
#' @param x 3-vectors or a 2-vectors (as rows of a matrix)
#'
#' Drops the 3 column if there is one
#' @export
drop3 <- function(x) {
  # truncate 3-vector to 2
  x[,-3,drop=FALSE]
}
#' Create a classless row vector (as a matrix)
#'
#' r is to row as c is to column
#' @export
r <- function(...) {
  # c is to column as r is to row
  rbind(c(...))
}
#' Normalize 3-vector
#'
#' @param x a 3-vector
#'
#' @export
n <- function(x) {
  if(ncol(x) != 3) stop('can only normalize 3-vectors')
  div <- ifelse(x[,3] == 0, 1, x[,3])
  x/div
}
#' Make 3-vector
#'
#' Only really makes sense for points
#'
#' @param x a 3-vector, 2-vector or vector of length 2 or 3.
#'
#' @export
to3 <- function(x, c3 = 1) {
  x <- rbind(x)
  if(ncol(x) == 2) cbind(x,c3)
  else if(ncol(x) != 3) stop('rbind(x) must have 2 or 3 columns')
  x
}
#' Make 2-vector
#'
#' Only really makes sense for finite points
#'
#' @param x a 3-vector, 2-vector or vector of length 2 or 3.
#'
#' @export
to2 <- function(x) {
  x <- rbind(x)
  if(ncol(x) == 2) return(x)
  else if(ncol(x) != 3) stop('rbind(x) must have 2 or 3 columns')
  x[,1:2,drop=FALSE]/x[,3]
}
#' Dual object
#'
#' @param o1,o2 two lines, or two points as 3-vectors
#'
#' @returns the intersection of two lines or the line joining two points
#'
#' @examples
#' dual(r(1,1,1), r(1,-1,1))
#'
#' @export
dual <- function(o1,o2) {
    ret <- pracma::cross(o1,o2)
    if('ppoint' %in% class(o1)) class(ret) <- 'pline'
    if('pline' %in% class(o1)) class(ret) <- 'ppoint'
    n(ret)
}
#' Points on a line
#'
#' @param x a 3-vector representing a line
#' @param len length of line, default 5
#'
#' Returns two points on a line as two rows of a matrix
#' as 2-vectors
#'
#' @export
topts <- function(line, len = 5) {
  line <- rbind(line)
  if(ncol(line) != 3) stop('line must be a 3-vector')
  pmid <- if(line[1,3] ==0) r(0,0,1) else dual(line,r(line[1,2],-line[1,1],0))
  pmid <- to2(pmid)
  orth <- r(line[1,2], - line[1,1])
  orth.len <- sqrt(sum(orth^2))
  p1 <- (pmid) + r(line[1,2],-line[1,1]) * len/orth.len
  p2 <- (pmid) - r(line[1,2],-line[1,1]) * len/orth.len
   rbind(p1,p2)
}
#' Translate a line
#'
#' See ![translating a line](drawings/translating_a_line.png)
#'
#' @param line represented as 3-vector
#' @param point represented as 3-vector or 2-vector
#'
#' @export
plus <- function(line,point) {
  point <- to2(point)
  ret <- line
  ret[,3] <- line[,3] - point[1] * line[,1] - point[2] * line[,2]
  ret
}
#' Intersection of line and 0-centered ellipse
#'
#' ![intersecting an ellipse](drawings/intersecting_an_ellipse.png)
#'
inter_ell_cline <- function(ell, line) {
  center <- attr(ell,'parms')$center %>% to2
  shape <- attr(ell,'parms')$shape
  radius <- attr(ell,'parms')$radius
  alphas <- chol(shape)%*%c(line)[1:2]
  theta <- atan2(-c(alphas)[1],c(alphas)[2])
  ret <- radius *t(chol(shape))%*%c(cos(theta),sin(theta))
  ret <- rbind(c(ret), -c(ret))
  ret
}
Sigma <- .1*diag(2)  + .8
ellobj <- ell(shape = Sigma)


inter_ell_cline(ellobj, r(1,1,0))

init <- function(from = -fac*c(1,1), to = fac*c(1,1), fac = 3,
                 xlab = '', ylab = '', asp = 1) {
  plot(rbind(from,to), type ='n', xlab = xlab, ylab = ylab, asp = asp)
}


init();abline(v=-10:10,col='gray50',lty=2);abline(h=-10:10,col='gray50',lty=2)

ellobj %>% lines
inter_ell_cline(ellobj, r(-1,2,0)) %>% lines
1/topts(r(1,1,1),10) %>% lines
topts(r(1,1,1) %>% plus(c(2,2)),10) %>% lines


topts(r(1,1,0),10) %>% lines
topts(r(1,0,2),10) %>% lines
topts(r(0,0,1),10) %>% lines


perp <- function(pline,through=r(0,0)) {
  # finds line that is perpendicular to pline
  # going through point through
  through <- as.ppoint(through)
  as.pline(pracma::cross(t(t(pline)*c(1,1,0)), through))
}

?ConjComp














ell_(c(1,1),shape = diag(2)+ .9 - .9*diag(2)) %>% plot(type= 'l', lwd = 2)
ell_(c(1,1),shape = diag(2)) %>% lines(col = 'red')
#
# Given two points draw a line
#
#' General:
#'
#' points are 1 x 2 matrices, created maybe with rbind
#' or 4 x 2 matrices i.e. intersection of two lines
#' lines are 2 x 2 matrices
#'
pline <- function(x,...) {
  UseMethod('pline')
}
pline.pline <- function(x,...) {
  x
}
pline.ppoint <- function(x1, x2, ...) {
  if(missing(x2)) {
    if(nrow(x1) < 2) stop('Need more than one row for x1 if no value for x2')
    x2 <- x1[-1,,drop=FALSE]
    x1 <- x1[rep(1,nrow(x1)-1),,drop = FALSE]
  }
  ret <- rbind(pracma::cross(x1,x2))
  class(ret) <- 'pline'
  ret
}
pline.default <- function(x,...) {
  # x must be an n x 2 matrix giving, usually, 2 points defining a line
  x <- rbind(x)
  if(ncol(x) != 2) stop('rbind(x) must have 2 columns')
  if(nrow(x) == 1) x <- rbind(c(0,0),x)
  ret <- cbind(x,1)
  class(ret) <- 'ppoint'
  pline(ret)
}
ppoint <- function(x,...) {
  UseMethod('ppoint')
}
ppoint.pline <- function(l1, l2,...) {
  if(missing(l2)) {
    if(nrow(l1) < 2) stop('Need more than one row for l1 if no value for l2')
    l2 <- l1[-1,,drop=FALSE]
    l1 <- l1[rep(1,nrow(l1)-1),,drop = FALSE]
  }
  ret <- rbind(pracma::cross(l1,l2))
  class(ret) <- 'ppoint'
  ret
}
as.pline <- function(x) {
  UseMethod('as.pline')
}
as.pline.pline <- function(x) {
  x
}
as.pline.default <- function(x) {
  x <- rbind(x)
  if(ncol(x) == 2) {
    x <- cbind(x,1)
  }
  if(ncol(x) != 3) stop('rbind(x) must have 2 or 3 columns')
  class(x) <- 'pline'
  x
}

as.ppoint <- function(x,...) {
  UseMethod('as.ppoint')
}

to2 <- function(x){
  x <- rbind(x)
  if(ncol(x) == 2) x
  else {
    x <- x/x[, 3]
    x[, -3, drop=FALSE]
  }
}
as.ppoint.default <- function(x) {
  x <- rbind(x)
  if(!(ncol(x) %in% c(2,3))) stop('rbind(x) must have 2 or 3 columns')
  if(ncol(x) == 2) x <- cbind(x,1)
  class(x) <- 'ppoint'
  x
}

perp(r(1,1,0))
perp(r(1,0,0))
perp(r(0,0,0))

perp(r(1,1,0), r(2,1))



init <- function(from = -fac*c(1,1), to = fac*c(1,1), fac = 3,
                 xlab = '', ylab = '', asp = 1) {
  plot(rbind(from,to), type ='n', xlab = xlab, ylab = ylab, asp = asp)
}
init()
# test

normalize <- function(x) {
  x <- rbind(x)
  if(ncol(x) != 3) stop('rbind(x) must have 3 columns')
  div <- ifelse(x[,3] == 0, 1, x[,3])
  x/div
}

to3 <- function(x,third) {
  x <- rbind(x)
  if(ncol(x) == 2) x <- cbind(x,third)
  x
}
to3(r(1,1), 1)
to3(r(1,1,1))




plot.pline <- function(x,...) {
  bounds <- par('usr')
  maxc <- max(abs(bounds))
  xx <- rbind(x)
  for(i in seq_len(nrow(xx))) {
    x <- xx[i,,drop=FALSE]
    disp(x)
    lp <- perp(x) # line perp to x thru 0
    disp(lp)
    ponline <- dual(x,lp)
    disp(ponline)
    ponline %>% to2 %>% points(pch=16,...)
    if(x[,3] != 0) { # line not through origin
      pmid <- ponline %>% to2
      p1 <- pmid + rev(pmid)*c(-1,1)
      p2 <- pmid - rev(pmid)*c(-1,1)
      p1m <- pmid + maxc * rev(pmid)*c(-1,1)
      p2m <- pmid - maxc * rev(pmid)*c(-1,1)
    } else {
      p1 <- x[1,2:1]*c(-1,1)
      p2 <- - x[1,2:1]*c(-1,1)
      p1m <- maxc * x[1,2:1]*c(-1,1)
      p2m <- - maxc * x[1,2:1]*c(-1,1)
    }
    rbind(p1,p2) %>% points(pch= 19,...)
    rbind(p1m, p2m) %>% lines(...)
  }
}

plot.ppoint <- function(x, ...) {
  points(to2(x),...)
}

plot.ppoint(rbind(r(1,1,1),r(3,3,2), r(3,3,0)))

r(3,3,0) %>% to2


init();abline(v=-10:10,col='gray50',lty=2);abline(h=-10:10,col='gray50',lty=2)
plot.pline(as.pline(r(.1,1,0)), col = 'red')
plot.pline(r(.1,1,1), col = 'blue')
plot.pline(as.pline(rbind(c(1,1,1),c(2,1,0), c(3,2,1))))
# plot(as.pline(rbind(c(1,1,1),c(2,1,0), c(3,2,1))))
plot.pline(rbind(c(1,1,1),c(2,1,0), c(3,2,1)))



plot.pline(as.pline(r(0,0,1)), col = 'blue')
plot.pline(as.pline(r(1,0,1)), col = 'blue')
plot.pline(as.pline(r(1,1,2)), col = 'blue')





plot.pline(as.pline(r(1,1,0)), col = 'red')
plot.pline(as.pline(r(1,1,1)), col = 'blue')


init()
plot.pline(as.pline(r(1,-1,0)), col = 'red')
plot.pline(as.pline(r(1,.1,0)), col = 'red')


abline(h=0)
abline(v=0)
perp(r(1,1,0))

plot.pline(as.pline(r(1,-1,1)), col = 'red')
plot.pline(perp(r(1,1,1)), col = 'blue')

perp(r(1,1,1))


rbindna <- function(x, ...) {
  if (nargs() == 0)
    return(NULL)
  if (nargs() == 1)
    return(x)
  rbind(x, NA, rbindna(...))
}

ellist <- function (center = rep(0, 2), shape = diag(2), radius = 1, n = 100,
                    angles = (0:n) * 2 * pi/n, fac = chol)
{
  rbindna <- function(x, ...) {
    if (nargs() == 0)
      return(NULL)
    if (nargs() == 1)
      return(x)
    rbind(x, NA, rbindna(...))
  }
  sx <- c(shape[1,1] %>% sqrt, 0)
  disp(sx)
  sy <- c(0, shape[2,2] %>% sqrt)
  disp(sy)
  yx <- c(sx[1], shape[1,2] / sqrt(shape[1,1]))
  disp(yx)
  xy <- c(shape[1,2] / sqrt(shape[2,2]), sy[2])
  ry <- c(0, (shape[2,2] - shape[1,2]^2/shape[1,1]) %>% sqrt)
  rx <- c((shape[1,1] - shape[1,2]^2/shape[2,2]) %>% sqrt, 0)
  co <- c(sx[1],sy[2])

  # if (missing(ellipse) && missing(diameters) && missing(box))
  #   all <- TRUE
  circle <- function(angle) cbind(cos(angle), sin(angle))
  rect <- cbind( c(-1,-1),c(-1,1),c(1,1), c(1,-1), c(-1,-1)) * sqrt(diag(shape))
  Tr <- fac(shape)
  ret <- list(ellipse = t(c(center) + t(radius * circle(angles) %*% Tr)),
              parallelogram = t(c(center) + t(radius * rbind(c(1, 1), c(-1, 1), c(-1, -1),
                                                             c(1,-1), c(1, 1)) %*% Tr)),
              box = t(center + (radius * rect)),
              yx = t(center + radius * cbind(yx, -yx)),
              xy = t(center + radius * cbind(xy, -xy)),
              diag1 = t(center + radius * cbind(co, -co)),
              diag2 = t(center + radius * cbind(co*c(1,-1), -co*c(1,-1))),
              ry = t(center + radius * cbind(ry, -ry)),
              rx = t(center + radius * cbind(rx, -rx))
  )

  # lapply(ret, dim) %>% print
  # do.call("rbindna", ret[c(ellipse, diameters, diameters, box, rectangle)])
  attr(ret,'ell') <- list(center = center, shape = shape)
  ret
}
#'
#' ![ell functions](ell.png)
#'
init()
ellist(shape = .1*diag(2) + .9) %>% lapply(lines)

setOldClass("foo")


setGeneric("type", function(x) standardGeneric("type"))
setMethod("type", signature("matrix"), function(x) "matrix")
setMethod("type", signature("character"), function(x) "character")
foo <- structure(list(x = 1), class = "foo")
type(foo)

?setOldClass
setGeneric('plus', function(e1,e2) standardGeneric('plus'))
setMethod("plus", signature(e1 = "foo", e2 = "ANY"), function(e1, e2) {
  disp('S4')
  structure(list(x = e1$x + e2), class = "foo")
})
setMetho
setMethod("type", signature("foo"), function(x) "foo")
showClass('foo')
type(foo)

plus(foo, 3)
foo %>% plus(3)
foo + 3

showMethods('+')

pracma::hypot

rms <- function(x) {
  M <- max(abs(x))
  if(M==0) return(0)
  x <- x/M
  M * sqrt(sum(x^2))
}
rms(1:100000)


library(pracma)
hy <- function(x,y) sqrt(x^2+y^2)
