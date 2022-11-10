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
  
  iter <- 1
  func_value <- func_value_original <- func(theta,...)
  grad_vector <- grad(theta,...)
  if (abs(func_value_original)==Inf) {
    stop('the objective function is not finite at the initial theta')
  }
  if (sum(abs(grad_vector))==Inf) {
    stop('the derivatives are not finite at the initial theta')
  }
  while (iter <= maxit) {
    grad_vector <- grad(theta,...)
    hess_matrix <- hess(theta,...)
    R <- chol(hess_matrix)
    inv_hess_matrix <- chol2inv(R)
    theta <- theta - drop(inv_hess_matrix %*% grad_vector)
    func_value <- func(theta,...)
    if (iter==max.half&func_value>=func_value_original) maxit <- maxit/2
    if (sum(abs(grad_vector)>=tol*abs(func_value+fscale))==0) {
      try(chol(hess_matrix))
      return(list(f=func_value,theta=theta,iter=iter,g=grad_vector,
                  Hi=inv_hess_matrix))
    } else {
      iter <- iter + 1
    }
  }
  warning('maxit is reached without convergence')
}