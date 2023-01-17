#'
#' THis is a part of a package for plotting objects
#' represented in projective 2 dimensional space
#' Using objects 'ppoint' and 'pline'
#'
# library(p3d)
library(spida2)
library(pracma)
?ellplus

drop3 <- function(x) {
  x[,-3,drop=FALSE]
}

drop3(r(1,1))

r <- function(...) {
  # c is to column as r is to row
  rbind(c(...))
}

ellplus(shape = diag(2)+ .5, box = T) %>% plot(type = 'l')
spida2::ellplus(shape = diag(2) + .5) %>% plot(type = 'l')
ellplus


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

perp <- function(pline,through=r(0,0)) {
  # finds line that is perpendicular to pline
  # going through point through
  through <- as.ppoint(through)
  as.pline(pracma::cross(t(t(pline)*c(1,1,0)), through))
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


dual <- function(o1,o2) {
  ret <- pracma::cross(o1,o2)
  if('ppoint' %in% class(o1)) class(ret) <- 'pline'
  if('pline' %in% class(o1)) class(ret) <- 'ppoint'
  ret
}
dual(as.ppoint(r(1,1,1)), r(1,-1,1)) %>% to2


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




