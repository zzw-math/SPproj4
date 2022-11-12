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
  #####
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
  #####
  func_value <- func(theta,...)   
  grad_vector <- grad(theta,...)
  if (abs(func_value)==Inf) {
    ## If the value of function approaches infinity, we stop the loop
    stop('the objective function is not finite at the initial theta')
  }
  if (sum(abs(grad_vector))==Inf) {
    ## If the gradient of function approaches infinity, we also stop the loop
    stop('the derivatives are not finite at the initial theta')
  }
  iter <- 0  ## initialize the iteration times
  while (iter <= maxit) {
    ## In while loop, we need to make sure for a given theta, the gradient
    ## matrix is positively definite. If it converges, we calculate Hessian to
    ## show positive definite matrix. If it diverges, we also calculate Hessian
    ## matrix, and then perform the iteration to get new theta.
    grad_vector <- grad(theta,...)
    if (!is.null(hess)) {
      hess_matrix <- hess(theta,...)
    } else {
      p <- length(theta)
      hess_matrix <- matrix(0,p,p)
      for (j in 1:p) {
        theta_add <- theta
        theta_add[j] <- theta_add[j]+eps
        hess_matrix[,j] <- (grad(theta_add,...)-grad_vector)/eps
      }
      hess_matrix <- (hess_matrix+t(hess_matrix))/2
    }
    if (max(abs(grad_vector))<tol*(abs(func_value)+fscale)) {
      ## In 'if', we test whether every gradient is lower than the critical
      ## value, and we define the number which greater than the critical value
      ## to be 0.
      ## When theta is satisfied, perform 'try' to test Hessian matrix.
      try(R <- chol(hess_matrix))
      inv_hess_matrix <- chol2inv(R)
      return(list(f=func_value,theta=theta,iter=iter,g=grad_vector,
                  Hi=inv_hess_matrix))
    }
    ## When theta is not satisfied, we replaced it by halving the step. 
    R <- chol(hess_matrix)
    elapse <- backsolve(R,forwardsolve(t(R),grad_vector)) 
    
    ## iteration elapse which will be halved
    theta_new <- theta - elapse 
    func_value_new <- func(theta_new,...)
    times.half <- 0  ## initialise the number of times to be halved
    while (func_value_new > func_value) {
      if (times.half==max.half) stop('the step fails to reduce the objective 
      despite trying max.half step halvings')
      ## Set the maximum of times we can halve the steps to 20
      elapse <- elapse/2  ## halved iteration elapse
      theta_new <- theta - elapse
      func_value_new <- func(theta_new,...)
      times.half <- times.half + 1 ## iterations of number of halves
    }
    iter <- iter + 1
    func_value <- func_value_new
    theta <- theta_new
  }
  warning('maxit is reached without convergence')
}

