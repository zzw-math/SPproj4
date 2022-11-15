##      names     | university user names
##   Ziwen Zhong  | s2326022
## Wenzheng Zhang | s2310185
##   Tianqi Dai   | s2302524

## github repo: https://github.com/zzw-math/SPproj4.git

## contribution
## Ziwen: 40% | Wenzheng: 30% | Tianqi: 30%
## Ziwen completed a draft code at first, Wenzheng and Tianqi added all the 
## comments. Then we tested the code and fixed a few errors together

#####
## summary
## Project 4 Overview: In this project, we define a function using Newton's
## Method to minimize the given function.
## Step 1: Check whether the function and gradient value is finite, if not, 
##         raise error.
## Step 2: Calculate gradient value and check if it is convergence. 
## Step 3: If it is convergence, judge the positive definite of Hessian, and 
##         return the results.
## Step 4: If it is not convergence, calculate the Hessian and step, then 
##         obtain new theta.
## Step 5: Check whether the step improves the object function.
## Step 6: If not, halve the step and back to step 4, we can only repeat this 
##         process max.half times or raise error.
## Step 7: Renew function value and theta and back to step 2. We can only repeat
##         maxit times or raise warning.
#####

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,
                 max.half=20,eps=1e-6) {
  #####
  ## Overview: function using Newton's Method to minimize the given function.  
  ## Input
  ##    theta: vector of float or integer. initial values of optimising 
  ##           parameters
  ##    func: function. the objective function to be minimised
  ##    grad: function. the function to calculate the gradient vector
  ##    hess: function. the function to calculate the Hessian matrix
  ##    tol: float. the convergence tolerance
  ##    fscale: integer. estimation of the magnitude of func near the optimum
  ##    maxit: integer. maximum number of Newton iterations
  ##    max.half: integer. maximum times to halve a step before the deduction 
  ##              of failure on objective improvement
  ##    eps: float. the finite difference intervals to use when a Hessian 
  ##         function is not provided 
  ## Output: Return error, warning or a list contains function value, theta, 
  ##         iteration times, gradient vector and inverse Hessian matrix (if 
  ##         have) at convergence.
  #####
  func_value <- func(theta,...)   
  grad_vector <- grad(theta,...)
  if (is.infinite(func_value)) {
    ## If the value of function approaches infinity, we stop the loop
    stop('the objective function is not finite at the initial theta')
  }
  if (any(is.infinite(grad_vector))) {
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
    ## If the function used to calculate the hessian matrix is provided, we can
    ## directly implement the hess function to get hessian matrix, otherwise,
    ## we need to estimate the hessian matrix by finite differencing of the 
    ## gradient vector.
    if (!is.null(hess)) {              
      hess_matrix <- hess(theta,...) 
    } else {
      p <- length(theta)
      hess_matrix <- matrix(0,p,p)   ## create an empty matrix to store hessian
      ## H_{.,j} = dg(theta)/d(theta_j), use for loop to generate each column.
      for (j in 1:p) {      
        theta_add <- theta 
        theta_add[j] <- theta_add[j]+eps ## theta_add add eps to j postition
        hess_matrix[,j] <- (grad(theta_add,...)-grad_vector)/eps
      }
      ## ensure the symmetry of hessian matrix
      hess_matrix <- (hess_matrix+t(hess_matrix))/2
    }
    ## Judge the convergence by the absolute of gradient value, whether they
    ## are all less than tol*(|func_value|+fsclae).
    if (max(abs(grad_vector))<tol*(abs(func_value)+fscale)) {
      ## Prepare a null object for inverse Hessian, try to calculate the
      ## Cholesky decomposition of Hessian, and then the inverse of Hessian.
      ## If failed to do this, which means Hessian is not positive definite,
      ## still return the result but the inverse Hessian is NULL.
      inv_hess_matrix <- NULL 
      try({  
        R <- chol(hess_matrix)    
        inv_hess_matrix <- chol2inv(R)  
      })     
      return(list(f=func_value,theta=theta,iter=iter,g=grad_vector,
                  Hi=inv_hess_matrix))
    }
    ## if not reach convergence, calculate the step (elapse) of iteration
    ## elapse = -H^{-1} %*% gradient
    R <- chol(hess_matrix)
    elapse <- -backsolve(R,forwardsolve(t(R),grad_vector)) 
    
    theta_new <- theta + elapse ## use step to obtain new theta 
    func_value_new <- func(theta_new,...) 
    times.half <- 0  ## initialize the number of times to be halved
    ## check whether the function value is lowered, if not, halve the step
    ## until it is lower or reaches max.half
    while (func_value_new > func_value) {
      if (times.half==max.half) stop('the step fails to reduce the objective 
      despite trying max.half step halvings')
      elapse <- elapse/2     ## halve the step
      theta_new <- theta + elapse
      func_value_new <- func(theta_new,...)
      times.half <- times.half + 1    ## halve times + 1
    }
    iter <- iter + 1    ## iteration times + 1
    func_value <- func_value_new    ## renew the function value
    theta <- theta_new              ## renew theta
  }
  ## If fail to reach convergence, raise a warning message.
  warning('maxit is reached without convergence') 
}

