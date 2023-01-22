#'
#' THis is a part of a package for plotting objects
#' represented in projective 2 dimensional space
#' Using objects 'ppoint' and 'pline'
#' 
#' General approach
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
#'
library(spida2)

n <- function(x) {
  if(ncol(x) != 3) stop('can only normalize 3-vectors')
  div <- ifelse(x[,3] == 0, 1, x[,3])
  x/div
}


pl <- function(x,...) {
  UseMethod('pl')
}

r <- function(...) {
  rbind(c(...))
}

l <- function(x,...) {
  UseMethod('l')
} 
l.line <- function(x,...) {
  x
}
l.point <- function(p1, p2,...) {
  if(missing(p2)) {
    if(nrow(p1) == 2) l(p(p1[1,]),p(p1[2,]))
    else l(l(p1,p(c(0,0,1))))
  } else {
    dual(p1,p2)
  }
}


l.default <- function(x,...) {
  x <- rbind(x,...)
  if(ncol(x) == 2) x <- cbind(x,1)
  x <- n(x)
  class(x) <- 'line'
  x
}

p <- function(x,...) {
  UseMethod('p')
} 
p.default <- function(x,...) {
  x <- rbind(x,...)
  if(ncol(x) == 2) x <- cbind(x,1)
  x <- n(x)
  class(x) <- 'point'
  x
}
p.line <- function(l1, l2, ...) {
  # intersection(s) of l1 and l2 
  if(missing(l2)) {
    if(nrow(l1) == 2) p(l(l1[1,]),l(l1[2,]))
    else stop('need 2 lines')
  } else {
    dual(l1,l2)
  }
}


pl.line <- function(x,...) {
  x <- l(x)
  for(i in seq_len(nrow(x))) {
    ll <- x[i,]
    if(ll[2] != 0) abline(a = -ll[3]/ll[2], b = -ll[1]/ll[2], ...)
    else abline(v = -ll[3]/ll[1], ...)
  }
  invisible(NULL)
}

pl.point <- function(x,...) {
  x <- p(x)
  disp(x)
  for(i in seq_len(nrow(x))) {
    ll <- x[i,,drop = FALSE]
    if(ll[3] == 1) points(ll[,1:2,drop=FALSE],...)
  }
  invisible(NULL)
}

pl.ell <-function(x, ...) {
  lines(x,...)
  invisible(NULL)
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



plus <- function(x,...) {
  UseMethod("plus")
}

plus.point <- function(obj,pt) {
  obj <- p(obj)
  pt <- p(pt)
  nr <- nrow(obj)
  ind <- seq_len(nrow(pt))
  ind <- rep(ind, length.out = nr)
  pt <- pt[ind,,drop=FALSE]
  ret <- obj * pt[,3] + pt * obj[,3]
  ret[,3] <- ret[,3]/2
  class(ret) <- 'point'
  ret
}

plus.line <- function(obj,pt) {
  obj <- l(obj)
  pt <- p(pt)
  nr <- nrow(obj)
  ind <- seq_len(nrow(pt))
  ind <- rep(ind, length.out = nr)
  pt <- pt[ind,,drop=FALSE]
  ret <- obj
  ret[,3] <- obj[,3] - pt[1] * obj[,1] - pt[2] * obj[,2]
  class(ret) <- 'line'
  ret
}

rows <- function(m) {
  ret <- lapply(seq_len(nrow(m)), function(i) {
    ret <- m[i,,drop=FALSE]
    class(ret) <- class(m)
    ret
  })
  ret
}

# rows(l(t(zm)))

dual <- function(x,...) {
  UseMethod('dual')
}
dual.line <- function(o1,o2) {
  # line through 2 points or point at intersection of two lines
  # - pipe through p or l if class uncertain
  a1 <- rbind(o1)
  a2 <- rbind(o2)
  # recycle o2 if fewer rows
  n1 <-nrow(a1)
  n2 <-nrow(a2)
  ind2 <- rep(seq_len(n2),length.out=n1)
  ret <- pracma::cross(a1,a2[ind2,])
  ret <- n(ret)
  if('point' %in% class(o1)) class(ret) <- 'line'
  if('line' %in% class(o1)) class(ret) <- 'point'
  # if(!missing(class.out)) class(ret) <- class.out
  ret
}
dual.point <- dual.line
dual.default <- dual.line
dual.ell <- function(x, ...) {
  shape <- attr(x,'params')$shape
  radius <- attr(x,'params')$radius
  ell(shape = solve(shape), radius = 1/radius)
}

ell <- function (center = rep(0, 2), shape = diag(2), radius = 1, n = 100,
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
  # disp(sx)
  sy <- c(0, shape[2,2] %>% sqrt)
  # disp(sy)
  yx <- c(sx[1], shape[1,2] / sqrt(shape[1,1]))
  # disp(yx)
  xy <- c(shape[1,2] / sqrt(shape[2,2]), sy[2])
  ry <- c(0, (shape[2,2] - shape[1,2]^2/shape[1,1]) %>% sqrt)
  rx <- c((shape[1,1] - shape[1,2]^2/shape[2,2]) %>% sqrt, 0)
  co <- c(sx[1],sy[2])
  
  # if (missing(ellipse) && missing(diameters) && missing(box))
  #   all <- TRUE
  circle <- function(angle) cbind(cos(angle), sin(angle))
  rect <- cbind( c(-1,-1),c(-1,1),c(1,1), c(1,-1), c(-1,-1)) * sqrt(diag(shape))
  Tr <- fac(shape)
  ret <- t(c(center) + t(radius * circle(angles) %*% Tr))
  params <- list(
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
  attr(ret,'params') <- c(params, list(center = center, shape = shape, radius = radius))
  class(ret) <- 'ell'
  ret
}

params <- function(obj) {
  attr(obj, 'params')
}


rot <- function(m,theta,...) {
  UseMethod('rot')
}
rot.default <- function(m,theta) {
  ret <- m %*% cbind(c(cos(theta), -sin(theta)), c(sin(theta),cos(theta)))
  attributes(ret) <- attributes(m)
  ret
}
rot.ell <- function(m,theta) {
  params <- attr(m,'params')
  rotm <- cbind(c(cos(theta), sin(theta)), c(-sin(theta),cos(theta)))
  shape <- rotm%*%params$shape%*%t(rotm)
  center <- params$center
  radius <- params$radius
  ell(center = center, shape = shape, radius = radius)
}
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

inter_ell <- (ell1,ell1) {
  # intersection of two centered ellipses
  par1 <- attr(ell1,'params')
  par2 <- attr(ell2,'params')
  cen <- sum(c(par1$center^2 ,par2$center^2))
  if(cen > 0) warning('intersections are for ellipses centered at origin')
  A <- solve(par1$radius^2 * par1$shape) - solve(par2$radius^2 * par2$shape)
  B <- solve(par1$radius^2 * par1$shape) 
  a <- solve(A) %*% B %*% cen
  k <- t(a)%*%A%*%a - t(cen)%*%B%*%cen
  


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




