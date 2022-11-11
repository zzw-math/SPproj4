source('proj4.r')


rb <- function(th,k=2) {
  k * (th[2]-th[1]^2)^2 + (1-th[1])^2
}

gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4 * th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

newt(theta=c(0,0),func=rb,grad=gb,hess=hb,k=2)


newt(theta=c(0,0),func=rb,grad=gb,k=2)


tol=1e-8;fscale=1;maxit=100;max.half=20;eps=1e-6
theta=c(0,0);func=rb;grad=gb
