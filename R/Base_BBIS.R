rbf_kernel <- function(x1,x2,h) {
  return( exp(-(1/h) * norm((x1-x2),type='2')**2) )
}

# x1 = px1 column vector
# x2 = px1 column vector
# s_x1 = px1 column vector of gradient at x1
# s_x2 = px1 column vector of gradient at x2
# h = shape parameter
# kernel = kernel function
b_rbf_kernel_p <- function(x1,x2,s_x1,s_x2,h,kernel) {
  k = kernel(x1,x2,h)

  kp1 = t(s_x1) %*% (k * s_x2)
  kp2 = t(s_x1) %*% (k * (2/h) * (x1-x2))
  kp3 = t(s_x2) %*% (k * (-2/h) * (x1-x2))
  kp4 = ((2*exp(-(x1-x2)**2)/h)%*%(h-2*(x1-x2)**2))/h**2

  return(kp1+kp2+kp3+kp4)
}


#' BBIS - Gradient Dependent Black Box Importance Sampling
#'
#' @param theta nxp matrix of sample points
#' @param theta_grads nxp matrix of sample point gradients
#' @param max_num_its integer of max number of optim steps
#' @param kernel type of kernel to use for importance sampling
#'
#' @return list of nx1 weights and 1xp weighted mean
#' @export
#'
#'
BBIS <- function(theta, theta_grad, max_num_its, kernel='rbf') {
  K_p = matrix(data=NA,nrow=nrow(theta),ncol=nrow(theta))

  dist = as.matrix(dist(theta, method = "euclidean",
                        diag = TRUE, upper = TRUE)**2)
  h = median(dist)

  for (ii in 1:nrow(K_p)) {
    for (jj in 1:ncol(K_p)) {
      if(kernel=='rbf') { K_p[ii,jj] = b_rbf_kernel_p(theta[ii,],theta[jj,],theta_grad[ii,],theta_grad[jj,],h,rbf_kernel) }
    }
  }

  weights = matrix(data=1,nrow=nrow(theta),ncol=1)

  f <- function(weights) {
    weights <- weights/sum(weights)

    t(weights)%*%K_p%*%weights
  }

  out = optim(weights,f,control=list(maxit = max_num_its,parscale=c(weights)),lower=0,method = 'L-BFGS-B')
  out$par = out$par/sum(out$par)


  ans = t(out$par) %*% theta
  new_points = matrix(data=NA,nrow=nrow(theta),ncol=ncol(theta))
  
  for (ii in 1:nrow(theta)) { new_points[ii,] = out$par[ii] * theta[ii,] }
  
  return(list(weights = out$par, adj_mean = ans, new_points = new_points))
}



