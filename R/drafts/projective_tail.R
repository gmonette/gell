


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


## axes ####

setGeneric('axes', function(x,y,type,tan,end,...) standardGeneric('axes'))
setOldClass('gell')
setOldClass('point')

## axes: gell point missing missing missing ####

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

setMethod("axes", signature(x = 'gell', y = 'point', type = 'missing', tan = 'logical',end = 'line'),
          function(x, y, type, tan, end){
            li <- axes(x, y, type = 'l', tan = tan)
            ret <- rbind(
              center = center(x),
              end = join(li, end))
            seg(p(ret))
          }
)


## axes: gell point character missing line ####

setMethod("axes", signature(x = 'gell', y = 'point', type = 'missing', tan = 'missing',end = 'line'),
          function(x, y, type, tan, end){
            axes(x, y, tan = FALSE, end = end)
          }
)

## axes: gell point character missing missing ####
setMethod("axes", signature(x = 'gell', y = 'point', type = 'character', tan = 'missing',end = 'missing'),
          function(x, y, type, tan){
            axes(x, y, type, FALSE)
          }
)



if(F){
  
  init(1.5)
  ge %>% pl
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

sel <- function(x, what, ...) {
  UseMethod('sel')
}

sel.list <- function(x, what,...) {
  x[[what]]
}

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

disp <- spida2::disp


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
init()
fun <- function(x,y) x^3 - y^2
inter_fun(fun,val = c(-1,-.1,0,.1,1), n =400) %>% pl()

pl.contour_list <- function(cl,...) {
  lapply(cl, function(x) {
    pl(x$contour, ...)
  }
  )
  invisible(cl)
}




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

setGeneric('inter', function(x, y, option,...) standardGeneric('inter'))

## inter: gell gell missing ####
setMethod("inter", signature(x = 'gell', y = 'gell', option = 'missing'),
          function(x, y, option, ...)  {
            inter_gells(x, y)
          }
)

## inter: gell gell numeric ####
setMethod("inter", signature(x = 'gell', y = 'gell', option = 'numeric'),
          function(x, y, option, ...)  {
            inter_gells(x, y, n=option)
          }
)

## inter: gell line missing ####
setMethod("inter", signature(x = 'gell', y = 'line', option = 'missing'),
          function(x, y, option, ...)  {
            inter_gell_line(x, y)
          }
)

## inter: gell seg missing ####
setMethod("inter", signature(x = 'gell', y = 'seg', option = 'missing'),
          function(x, y, option, ...)  {
            inter_gell_seg(x, y)
          }
)

## inter: gell point missing ####
setMethod("inter", signature(x = 'gell', y = 'point', option = 'missing'),
          function(x, y, option, ...)  {
            inter_gell_seg(x, p(y))
          }
)


rmat <- function(rho) {
  cbind(c(1,rho), c(rho, 1))
}


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

setGeneric('fta', function(x,y,z,...) standardGeneric('fta'))

## fta point point missing ----------
setMethod("fta", signature(x = 'point', y = 'point', z = 'missing'),
          function(x, y, z, ...)  {
            seg(rbind(x,y))
          }
)

## fta point line missing ----------
setMethod("fta", signature(x = 'point', y = 'line', z = 'missing'),
          function(x, y, z, ...)  {
            to <- join(perp(y,x),y)
            seg(rbind(x,to))
          }
)

## fta point point missing ----------
setMethod("fta", signature(x = 'line', y = 'point', z = 'missing'),
          function(x, y, z, ...)  {
            from <- join(perp(x,y),x)
            seg(rbind(from,y))
          }
)

## fta line line line ----------
setMethod("fta", signature(x = 'line', y = 'line', z = 'line'),
          function(x, y, z, ...)  {
            # z is used as a direction
            from <- join(x,z)
            to <- join(y,z)
            seg(rbind(from,to))
          }
)

## fta line line point ------------------------
setMethod("fta", signature(x = 'line', y = 'line', z = 'point'),
          function(x, y, z, ...)  {
            stop('no method from line to line along point')
          }
)

## fta line point line ------------------------------
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
setMethod("fta", signature(x = 'point', y = 'line', z = 'line'),
          function(x, y, z, ...)  {
            ret <- fta(y, x, z)
            seg(ret[2:1,])
          }
)


## fta line point point ---------------------------
setMethod("fta", signature(x = 'line', y = 'point', z = 'point'),
          function(x, y, z, ...)  {
            # z used as a direction
            z <- join(z,p(0,0,1))
            fta(x,y,z)
          }
)

## fta point line point  -------------------------------
setMethod("fta", signature(x = 'point', y = 'line', z = 'point'),
          function(x, y, z, ...)  {
            z <- join(z,p(0,0,1))
            fta(x,y,z)
          }
)

init()
gell(shape = rmat(.8), center = c(1,1)) %>% pl -> ge
ge %>% axes(p(0,1,0)) %>% pl
ge %>% axes(p(0,1,0)) %>% sel(2) %>%  join(p(0,0,1)) %>% pl
ge %>% axes(p(0,1,0)) %>% {join(.$points,.$tan_dirs)} %>% pl
ge %>% axes(p(0,1,0))



## fta gell line point -------------------------------------
setMethod("fta", signature(x = 'gell', y = 'line', z = 'point'),
          # use point a direction towards line
          # so only the center matters
          function(x, y, z, ...)  {
            from <- p(x$center)
            fta(from,y,z)
          }
)

## fta gell line line -------------------------------------
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
setMethod("fta", signature(x = 'gell', y = 'point', z = 'missing'),
          function(x, y, z, ...)  {
            fta(p(x$center),y)
          }
)

comb <-function(g1,g2,lam = .5) {
  w1 <- (1-lam)  * g1$gamma %*% diag(1/g1$d^2) %*% t(g1$gamma)/g1$radius
  w2 <- (lam) * g2$gamma %*% diag(1/g2$d^2) %*% t(g2$gamma)/g2$radius
  cen <- w1 %*% g1$center + w2 %*% g2$center
  var <- solve(w1 + w2)
  cen <- var %*% cen
  gell(shape = var, center = cen)
}

## fta gell gell missing
setMethod("fta", signature(x = 'gell', y = 'gell', z = 'missing'),
          function(x, y, z, ...)  {
            ret <- sapply(seq(0,1, length.out = 1001) , 
                          function(lam) comb(x, y, lam)$center)
            seg(p(t(ret)))
          }
)
setMethod("fta", signature(x = 'gell', y = 'gell', z = 'numeric'),
          function(x, y, z, ...)  {
            ret <- sapply(seq(0,1, length.out = z) , 
                          function(lam) comb(x, y, lam)$center)
            seg(p(t(ret)))
          }
)

if(F){
  
  fta(g1,g2,5)  %>% l %>% pl(pch = 16,cex = 2 , col = 'red')
  fta(center(g1),l(0,1,-4))  %>% l %>% pl(pch = 16,cex = 2 , col = 'red')
  l(0,-1,4) %>% fta(fta(g1,g2,5) %>% p,.) %>% pl
}


## ptr generic for projective transformations  -------------------------
setGeneric('ptr', function(x,y,z,...) standardGeneric('ptr'))

## ptr point matrix missing ----------
setMethod("ptr", signature(x = 'point', y = 'matrix', z = 'missing'),
          function(x, y, z, ...)  {
            p(t(y %*% t(x)))
          }
)

## ptr line matrix missing ----------
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




if(F){
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
  ge %>% axes(p(1,2,0)) 
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
  gell.default
}

COND <- FALSE
if(COND) {
  init()
  
  # draw an ellipse
  
  init()
  abline(h=0)
  abline(v=0)
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
    init()
    gell(shape = rmat(.8)) %>% pl -> g1
    gell(shape = rmat(-.7), center = c(4,-1)) %>% pl -> g2
    
    
    comb(g1,g2) %>% pl
    fta(g1,g2) %>%  pl
    
    # How evidence overwhelms a prior
    
    init()
    gell(radius=3) %>% pl(col = 'red')->gprior
    p(2,-2,1) %>% text(lab('prior'), col = 'red')
    
    freq_ell <- lapply(c(2,5,10,25,50,100),
                       function(n) {
                         gell(shape=rmat(-.9), center = c(3,4), radius = 1/sqrt(n))
                       }
    )
    freq_ell %>% pl
    post_ell <- lapply(freq_ell, function(ge) comb(ge,gprior))
    post_ell %>% pl(col='green')
    
    
    
    
    
    for(lam in seq(0,1,.01)) {comb(g1,g2,lam) %>% {c(.$center)} %>%  p %>% pl}
    for(lam in seq(0,1,.01)) {comb(g1,g2,lam) %>%  pl(col = '#77777722')}
    comb(g1,g2)    %>% pl (fill = lightgray)
    fta(g1,g2) %>%  pl(col = hcl(0,50,100,.5))
  }
  
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


x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
cl <- contourLines(x, y, volcano)
## summarize the sizes of each the contour lines :
cbind(lev = vapply(cl, `[[`, .5, "level"),
      n  = vapply(cl, function(l) length(l$x), 1))

z <- outer(-9:25, -9:25)
pretty(range(z), 10) # -300 -200 ... 600 700
utils::str(c2 <- contourLines(z))
# no segments for {-300, 700};
#  2 segments for {-200, -100, 0}
#  1 segment  for  100:600

?contourLines




