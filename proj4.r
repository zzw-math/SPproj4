##      names     | university user names
##   Ziwen Zhong  | s2326022
## Wenzheng Zhang |
##   Tianqi Dai   | s2302524

## github repo: https://github.com/zzw-math/SPproj4.git

## contribution
## 

#####
## summary
## Project 4 Overview: In this project, we are going to define a function newt.
#####

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,
                 max.half=20,eps=1e-6) {
  ## function newt with default data  
  ## Input: theta, vector of the initial values of optimising parameters
  ##        func, the objective function to be minimised
  ##        grad, the gradient vector returned by the gradient function
  ##        hess,  the Hessian matrix function
  ##        tol, the convergence tolerance
  ##        fscale, estimation of the magnitude of func near the optimum
  ##        maxit, maximum number of Newton iterations
  ##        max.half,maximum times to halve a step before the deduction of
  ##                failure on objective improvement
  ##        eps, the finite difference intervals to use when a Hessian function
  ##                is not provided 
  
  iter <- 0
  func_value <- func(theta,...)
  grad_vector <- grad(theta,...)
  if (abs(func_value)==Inf) {
    stop('the objective function is not finite at the initial theta')
  }
  if (sum(abs(grad_vector))==Inf) {
    stop('the derivatives are not finite at the initial theta')
  }
  while (iter <= maxit) {
    times.half <- 0
    grad_vector <- grad(theta,...)
    if (sum(abs(grad_vector)>=tol*abs(func_value+fscale))==0) {
      hess_matrix <- hess(theta_new,...)
      try(chol(hess_matrix))
      R <- chol(hess_matrix)
      inv_hess_matrix <- chol2inv(R)
      return(list(f=func_value,theta=theta_new,iter=iter,g=grad_vector,
                  Hi=inv_hess_matrix))
    } else {
      hess_matrix <- hess(theta,...)
      R <- chol(hess_matrix)
      elapse <- backsolve(R,forwardsolve(t(R),grad_vector))
      theta_new <- theta - elapse
      func_value_new <- func(theta,...)
      while (func_value_new>func_value) {
        if (times.half==20) stop('the step fails to reduce the objective despite 
                               trying max.half step halvings')
        elapse <- elapse/2
        theta_new <- theta - elapse
        func_value_new <- func(theta_new,...)
        times.half <- times.half + 1
      }
      iter <- iter + 1
      func_value <- func_value_new
      theta <- theta_new
    }
  }
  warning('maxit is reached without convergence')
}
