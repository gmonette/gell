# Testing stuff
# 
# order of points on a line
# 
p(2,3) %>% rbind(p(1,-1)) %>% p

p(1,1,1) %>% plus(p(2,1,1))
p(1,1,1) + p(2,1,1)


p(2,3) %>% join(p(1,-1)) 
showMethods(join)

xp <- function(x,y) {
  pracma::cross(x,y)
}

b1 <- p(1,1,1)
b2 <- p(3,2,1)
xp(b1,b2)
a <- xp(b1, xp(b2,b1))
a
dp <- function(x,y) {
  sum(x*y)
}

dp(a, b1)
dp(a, b2)

wh <- function(b1,b2) {
  b1 <- p(b1)
  b2 <- p(b2)
  a <- xp(b1, xp(b2,b1))
  a <- a/dp(a,b2)
  a
}

wh <- function(x, b1, b2, intercept = F) {
  b1 <- p(b1)
  b2 <- p(b2)
  x <- p(x)
  f <- lsfit(t(rbind(b1,b2)), t(x), intercept = F)
  f
}

b1 <- p(1,1,1)
b2 <- p(2,3,1)
x <- p(0,-1,1)
X <- rbind(
  p(1.5,2,1),
  p(2,3,1),
  p(0,-1,1),
  p(1,2,0)
)

wh(p(1.5, 2, 1), b1, b2)
wh(X, b1, b2)

wh(p(1.5, 0, 1), p(1,0,1), p(2,0,1))$coef
wh(p(12, 0, 1), p(1,0,1), p(2,0,0))$coef

wh(p(1.5, 0, 1), p(1,0,1), p(2,0,1))$coef

wh(p(1.5, 0, 1), p(1,0,1), p(2,0,1))$coef

wh(p(1.5, 0, 1), p(1,0,1), p(2,0,1))$coef
wh(p(1.5, 0, 1), p(1,0,1), p(2,0,0))[c('coefficients','residuals')]

wh(p(1000, 0, 1), p(1,0,1), p(3,0,1))[c('coefficients','residuals')]

X <- cbind(seq(-1,11), 0, 1)
rownames(X) <- paste0('X:',seq(-1,11))

wh(X, p(1,0,1), p(10,0,1))[c('coefficients','residuals')]
wh(X, p(1,0,1), p(10,0,0))[c('coefficients','residuals')]



wh(b1,b2) %>% dp(b1)
wh(b1,b2) %>% dp(b2)
b1
b2
bn <- b1 + 10 * (b2 - b1)
wh(b1,b2) %>% dp(p(bn))




wh(b1,b2) %>% dp(p(bn))

wh(b1,b2) %>% dp(p(bn))
wh(b1,b2) %>% dp(b2-b1)
wh(b1,b2) %>% dp(b1 -b2)
bn <- b1 - 10 * (b2 - b1)
wh(b1,b2) %>% dp(bn)

b1
b2 <- p(2,3,0)

wh(b1,b2)
wh(b1,b2) %>% dp(b1)
wh(b1,b2) %>% dp(b2)
wh(b1,b2) %>% dp(b1 + 10*b2)

wh(b1,b2) %>% dp(b1-10*b2)




