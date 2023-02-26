#'
#' This is a part of a package for plotting objects
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
## FUNCTIONS ####

library(pracma)
library(spida2)



disp <- spida2::disp

init <- function(from = -fac*c(1,1), to = fac*c(1,1), fac = 3,
                 xlab = '', ylab = '', asp = 1) {
  plot(rbind(from,to), type ='n', xlab = xlab, ylab = ylab, asp = asp)
}

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


COLS <- c('red','cyan','purple','yellow','teal' ,'blue','#4444FF','green','goldenrod',
               '#23423411',
               '#014488')
red <- 'red'
  blue <- 'blue'
    green <- 'green'
      black <- 'black'
               

#' Normalize homogeneous coordinates
#' 
#' @param x an object represented as a n x 3 numerical matrix
#' 
#' @export
n <- function(x) {
  if(ncol(x) != 3) stop('can only normalize homogeneous 3-vectors')
  div <- ifelse(x[,3] == 0, pracma::hypot(x[,1],x[,2]), x[,3])
  ret <- x/div
  class(ret) <- class(x)
  ret
}
#' Plot generic function

pl <- function(x, ...) {
  # plot method
  UseMethod('pl')
}

r <- function(...) {  # maybe get rid of this
  rbind(c(...))
}

##+ line ####
setOldClass('line')
#' Constructor for line object
#' 
#' Lines are represented by rows of n x 3 matrices where
#' each row is the projective representation of a 2D line
#' in homogeneous coordinates.
#' 
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
#' @param x an object 
#' 
#' @returns
#' A line object consisting of an n by 3 matrix where each row
#' represents a line.
#' 
#' @export
l <- function(x,...) {
  # constructor for 'line'
  UseMethod('l')
} 


l.line <- function(x,...) {
  n(x)
}

l.default <- function(x, ...) {
  if(length(x) == 1) x <- rbind(c(x, ...))  # why???
  else x <- rbind(x,...)
  if(ncol(x) == 2) x <- cbind(x,0)  # line through origin
  x <- n(x)           # normalize
  class(x) <- 'line'
  x
}


##+ point constructors ####

setOldClass('point')


p <- function(x, ...) {
  UseMethod('p')
} 

p.default <- function(x, ...) {
  if(length(x) == 1) x <- rbind(c(x, ...))
  else x <- rbind(x,...)
  if(ncol(x) == 2) x <- cbind(x,1)
  x <- n(x)
  class(x) <- 'point'
  x
}

##' join #### 

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

# test
join(p(c(1,1),c(1,-1)), p(c(2,2)))
join(
  p(c(1,1,0)),  # point at infinity)
  p(c(1,2))
)
join(p(c(1,1,1)), p(c(1,1,1))) %>% join(l(c(1,-1,0))) %>% join(p(c(1,1,1)))

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

##' Plotting methods for lines and points


pl <- function(x,...) {
  UseMethod('pl')
}

pl.default <- function(x, ...) {
  cat(paste("\\npl: no method for class:", class(x), "\\n"))
  invisible(x)
}

pl.matrix <- function(x, type = 'p', ...) {
  points(x[,1:2, drop = FALSE], type = type, ...)
  invisible(x)
}

pl.line <- function(x,...) {
  x <- l(x)
  for(i in seq_len(nrow(x))) {
    ll <- x[i,]
    if(ll[2] != 0) abline(a = -ll[3]/ll[2], b = -ll[1]/ll[2], ...)
    else abline(v = -ll[3]/ll[1], ...)
  }
  invisible(x)
}

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

# test
init <- function(fac= 5) plot(fac*c(-1,1),fac*c(-1,1), type = 'n', asp = 1)
if(F){
  
init()
p(c(0,0)) %>% 
  join(p(c(1,1))) %>% pl -> line1
p(c(1,1)) %>% 
  join(p(c(-1,1))) %>% pl -> line2
line1 %>% join(line2) %>% pl(pch=16, col = 'red')

}

## plus, mult ####
 
setGeneric('plus', function(x,y,...) standardGeneric('plus') )
setGeneric('mult', function(x, y, ...) standardGeneric('mult') )

setMethod('mult', signature('line','numeric'),
          function(x, y, ...) {
            x <- l(x)
            l(cbind(x[,1:2,drop=FALSE]/y,x[,3]))
          } 
)

setMethod('mult', signature('line','matrix'),
          function(x, y, ...) {
            # matrix operates on the left on points in line as a column vector 
            x <- l(x)
            A <- solve(y)
            if(nrow(y) == 2 && ncol(y) ==2) A <- cbind(rbind(A, 0), c(0,0,1))
            l(x %*% A)
          } 
)

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

rotm <- function(theta = pi/2) {
  cbind( c(cos(theta), sin(theta)), c(-sin(theta), cos(theta)))
}


rotm(pi/100)
##### NEXT: circle intersects line ####
##### then segments                ####



if(F){
  
  init()
  p(0,0,1) %>% pl
  l(1,1,1) %>% pl
  l(1,1,1) %>% mult(rotm(.1)) %>% pl
  l(1,1,1) %>% mult(rotm(pi/3)) %>% pl
  p(2,1,1) %>% {mult(rotm(pi/3),.)} %>% pl
  p(2,1,1) %>% {mult(rotm(0/3),.)} %>% pl
  
  
  l(1,1,0) %>% pl
  l(1,1,0) %>% mult(cbind(c(1,0),c(2,1))) %>% pl
  l(1,1,0) %>% mult(cbind(c(1,0),c(3,1))) %>% pl
  p(0,0,1) %>% pl(pch = 16)
  l(1,1,1) %>% pl
  l(1,1,1) %>% mult(cbind(c(1,0),c(3,1))) %>% pl
  l(1,1,1) %>% mult(cbind(c(0,-1),c(1,0))) %>% pl
}


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

setMethod('plus', signature('point','line'),
  function(x, y, ...) {
    plus(y, x, ...)
  } 
)

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





if(FALSE){
  
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
}

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
  # turn each row into element of list
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


circle <- function(n = 100) {
  theta <- seq(0, 2*pi, length.out = n)
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

params <- function(obj) {    # not needed ??????
  attr(obj, 'params')
}

##' rot ####
##' 
rot <- function(m,theta,...) {
  UseMethod('rot')
}
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

rot.point <- function(m,theta) {
  ret <- m %*% cbind(c(cos(theta), -sin(theta),0), c(sin(theta),cos(theta),0),c(0,0,1))
  attributes(ret) <- attributes(m)
  ret
}
rot.line <- function(m,theta) {
  ret <- m %*% cbind(c(cos(theta), -sin(theta),0), c(sin(theta),cos(theta),0),c(0,0,1))
  attributes(ret) <- attributes(m)
  ret
}
rot.seg <- function(m,theta) {
  ret <- m %*% cbind(c(cos(theta), -sin(theta),0), c(sin(theta),cos(theta),0),c(0,0,1))
  attributes(ret) <- attributes(m)
  ret
}

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
setOldClass('gell')
#' 
#' Constructor for gell objects
gell <- function(x,...) {
  if(missing(x)) gell.default(...)
  else UseMethod('gell')
}

# gell.ell <- function(ell,...) {
#   params <- attr(m,'params')
#   ei <- eigen(params$shape)
#   ret <- list(gamma = ei$vectors, d = sqrt(ei$values), radius = params$radius, center = params$center)
#   class(ret) <- 'gell'
#   ret
# }

gell.default <- function(
    center = c(0,0), shape = diag(2), radius = 1, 
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
if(FALSE) {
init()  
  gell(A=cbind(c(2,1),c(0,4))) %>%  pl(col = 'red')
  gell(A=cbind(c(2,1),c(0,4))) %>% rot(pi/4) %>% pl(col = 'red')
}


gell.matrix <- function(X,...) {
  # data ellipse
  center <- apply(as.matrix(X), 2, mean)
  X <- t(t(X) - center)/sqrt(nrow(X))
  sv <- svd(X, nu = 0)
  ret <- list(gamma = sv$v, d = sv$d, radius = 1, center = center )
  class(ret) <- 'gell'
  ret
}  

pl.gell <- function(gell, ...) {
  if(max(gell$d) < Inf & min(gell$d) > 0) {
    U <- circle()
    pts <- t(gell$gamma %*% (t(U) * gell$radius * gell$d ) + gell$center)
    lines(pts, ...)
  }
  else {
    warning('pl.gell: case not yet implemented')
  }
  invisible(gell)
}

center <- function(gell,...) {
  # generic in spida2
  p(gell$center)
}
pcenter <- center

if(FALSE) {
  
  X <- cbind(rnorm(200),rnorm(200)) %*% cbind(c(1,1), c(.2,-.3))
  init()
  X %>% p %>% pl  
  gell(X) %>% pl
  methods(pl)
  
  gell(X) %>% center %>% pl(pch = 16, col = 'red')
  
}

## axes ####

setGeneric('axes', function(x,y,...) standardGeneric('axes'))
setOldClass('gell')
setOldClass('point')

setMethod("axes", signature(x = 'gell', y = 'point'),
  function(x, y, ...)  {
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

sel <- function(x, nam) x[[nam]]

setOldClass('gell')
setOldClass('point')
setOldClass('line')

if(FALSE) {
  
  X <- matrix(rnorm(20), ncol = 2)
  X <- t(c(1,2)+t(X))
  init()
  X %>% pl
  gell(X) %>% pl
  gell(X) %>% axes(p(1,1)) %>% sel('points') %>% pl(col = 'red')
  gell(X) %>% axes(p(1,1)) %>% sel('tan_dirs') %>% join(center(gell(X))) %>% pl
  gell(X) %>% center %>% pl(col = 'red')
  gell(X) %>% axes(p(1,0,0)) %>% pl
}

# gell_conj <- function(ge, dir) {
#   # line thru center conjugate to a direction
#   # i.e. intersecting ellipse with tangent parallel to a direction
#   a <- axes(ge, dir)
# }

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

if(F){
  init()
  l(1,1,1) %>% pl
  l(1,1,1) %>% plus(p(1,1,1))%>% pl
  p(0,0,1) %>% pl(pch = 16)
  p(1,1,1) %>% pl(pch = 16)
  
  gell(X) %>% pl
  gell(X) %>% gell_box(p(1,1,0)) %>% pl
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

pl.list <- function(x, ...) {
  lapply(x, pl, ...)
  invisible(x)
}

#######   Intersections     #####################


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

inter_gell_seg <- function(ge, sg) {
  # %>% HIGHLY %>% BROKEN
 # maybe we son't need this is including circle includes all finite points
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






if(F){
  

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
#############  SEGMENTS  #################################

##' DO LATER ####

# pl.ell <-function(x, ...) {
#   lines(x,...)
#   invisible(x)
# }

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

seg.default <- seg.point

close.seg <- function(sg) {
  if(!identical(sg[1,],sg[nrow(sg),])) sg <- rbind(sg, sg[1,])
  seg(p(sg))
}

seg.point <- function(pp, close = FALSE, ...) {
  pp <- p(pp)
  if(close) pp <- close.seg(pp)
  class(pp) <- 'seg'
  pp
}

seg.default <- seg.point

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

if(F){
  
l(rbind(c(1,1,1), l(1,-1,1)))
join(l(1,1,1), l(1,-1,1)) %>% 
  rbind(p(2,3,1)) %>% p %>%  seg
}

pl.seg <- function(x, close = FALSE, type = 'l', fill = NULL, ...) {
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
  disp(rad)
  
  # ge <- gell(radius = rad)  # for points at infinity
  
  if(close) x <- close.seg(x)
  
  
  
  if(x[1,3] == 0) {
    ret <- matrix(0, ncol = 3, nrow =0)
  } else {
    ret <- x[1,,drop = FALSE]
  }
  
  att <- list()     
  for(i in 2:nrow(x)) {
    disp(i)
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
        disp('fin-fin')
        ret <- rbind(ret, x[i,])
        disp(ret)
      } else {                             # previous point infinite
        disp('infin-fin')
        ge <- gell(center=x[i,1:2], radius = rad)
        disp(ge)
        att <- append(att,list(ge))
        infp <- p(axes(ge, p(x[i-1,]))[[1]][1,])
        att <- append(att,list(infp))
        ret <- rbind(ret, infp, ppi = x[i,])
        disp(ret)
      }
    } else {  # current point infinite
      
      if(x[i-1,3] == 0) {  # previous point infinite
      disp('infin-infin')
      } else {  # previous point finite
        disp('fin-infin')
        
        ge <- gell(center=x[i-1,1:2], radius = rad)
        disp(ge)
        att <- append(att, list(ge))
        infp <- p(axes(ge, p(x[i,]))[[1]][1,])
        att <- append(att, list(infp))
        ret <- rbind(ret, infp)
        disp(ret)

      }
    }
  }
  if(is.null(fill)) {
     lines(ret[,1:2, drop = FALSE], type = type, ...)
  } else {
     polygon(ret[,1],ret[,2], col = fill, ...)
  }
  disp(att)  
  ret <- seg(p(ret))
  attr(ret,'att') <- att
  invisible(ret)
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
               
        




