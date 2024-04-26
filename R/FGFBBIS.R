rbf_kernel <- function(x1,x2,h) {
  return( exp(-(1/h) * norm((x1-x2),type='2')) )
}


# x1 and x2 are points
# px2 is the probability of x2, p(x2)
# kernel is the rbf kernel function
f_gf_rbf_kernel_p <- function(x1,x2,px1,px2,h,kernel) {
  k = rbf_kernel(x1,x2,h)

  kp1 = sum(k * (2/h) * (x1-x2))
  kp2 = sum(k * (-2/h) * (x1-x2))
  kp3 = (k * (2/h)) + ((-2/h) * (x1-x2) %*% (k * (2/h)* (x1-x2)))

  return((k+kp1+kp2+kp3)*(px1*px2))
}


#' FGFBBIS - Full Gradient Free Black Box Importance Sampling
#'
#' @param theta nxp matrix of sample points
#' @param theta_probs nx1 vector of point probabilities
#' @param max_num_its integer of max number of optim steps
#' @param kernel type of kernel to use for importance sampling
#'
#' @return list of nx1 weights and 1xp weighted mean
#' @export
#'
#'
FGFBBIS <- function(theta, theta_probs, max_num_its, kernel='rbf') {
  K_p = matrix(data=NA,nrow=nrow(theta),ncol=nrow(theta))

  dist = as.matrix(dist(theta, method = "euclidean",
                        diag = TRUE, upper = TRUE))
  h = median(dist)

  for (ii in 1:nrow(K_p)) {
    for (jj in 1:ncol(K_p)) {
      if(kernel=='rbf') { K_p[ii,jj] = f_gf_rbf_kernel_p(theta[ii,],theta[jj,],theta_probs[ii,],theta_probs[jj,],h,rbf_kernel) }
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

  return(list(weights = out$par, adj_mean = ans))
}

