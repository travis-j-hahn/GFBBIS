rbf_kernel <- function(x1,x2,h) {
  return( exp(-(1/h) * norm((x1-x2),type='2')) )
}


# x1 = px1 column vector
# x2 = px1 column vector
# s_x1 = px1 column vector of gradient at x1
# s_x2 = px1 column vector of gradient at x2
# px1 = probability of x1
# px2 = probability of x2
# sub_px1 = probability of x1 from surrogate distribution
# sub_px2 = probability of x2 from surrogate distribution
# h = shape parameter
# kernel = kernel function
gf_rbf_kernel_p <- function(x1,x2,s_x1,s_x2,px1,px2,sub_px1,sub_px2,h,kernel) {
  k = kernel(x1,x2,h)

  # print(sub_px1*sub_px2)
  # print(px1*px2)

  kp1 = t(s_x1) %*% (k * s_x2)
  kp2 = t(s_x1) %*% (k * (2/h) * (x1-x2))
  kp3 = t(s_x2) %*% (k * (-2/h) * (x1-x2))
  kp4 = (k * (2/h)) + ((-2/h) * (x1-x2) %*% (k * (2/h)* (x1-x2)))

  return((kp1+kp2+kp3+kp4)*((sub_px1*sub_px2)/(px1*px2)))
}


# Returns multivariate normal log probabilities and grad log probs
# theta = nxp matrix of points
# mu = px1 vector of mean
# sigma = pxp covariance matrix
mvtnorm_prob <- function(theta,mu,sigma) {
  probs = matrix(data=NA,nrow=nrow(theta),ncol=1)
  grad_probs = matrix(data=NA,nrow=nrow(theta),ncol=ncol(theta))


  for (ii in 1:nrow(theta)) {
    probs[ii] = ((2*pi)**(-ncol(theta)/2)) * (det(sigma)**(-1/2)) * exp((-1/2)*(t(theta[ii,]-mu)%*%solve(sigma)%*%(theta[ii,]-mu)))

    grad_probs[ii,] = -solve(sigma)%*%(theta[ii,]-mu)
  }

  return(list(den=(probs),grad_x=grad_probs))
}



#' GFBBIS - Gradient Free Black Box Importance Sampling, surrogate density is mvt normal
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
GFBBIS <- function(theta, theta_prob, max_num_its, kernel='rbf') {
  # create new multi-var distribution from data and get gradients from that
  mvtnormdata = mvtnorm_prob(theta,mu=colMeans(theta),sigma=cov(theta)*3)
  sub_theta_prob = mvtnormdata$den
  sub_theta_grad = mvtnormdata$grad_x

  print(sum(sub_theta_prob))

  dist = as.matrix(dist(theta, method = "euclidean",
                        diag = TRUE, upper = TRUE)**2)
  h = median(dist)

  K_p = matrix(data=NA,nrow=nrow(theta),ncol=nrow(theta))

  for (ii in 1:nrow(K_p)) {
    for (jj in 1:ncol(K_p)) {
      if(kernel=='rbf') { K_p[ii,jj] = gf_rbf_kernel_p(theta[ii,],theta[jj,],sub_theta_grad[ii,],sub_theta_grad[jj,],
                                                      theta_prob[ii,],theta_prob[jj,],sub_theta_prob[ii,],sub_theta_prob[jj,],
                                                      h,rbf_kernel) }
    }
  }

  weights = matrix(data=1,nrow=nrow(theta),ncol=1)

  f <- function(weights) {
    weights <- weights/sum(weights)

    t(weights)%*%K_p%*%weights
  }

  out = optim(weights,f,control=list(maxit = max_num_its,parscale=c(weights)),lower=0,method = 'L-BFGS-B')
  # print(out)
  out$par = out$par/sum(out$par)


  ans = t(out$par) %*% theta

  return(list(weights = out$par, adj_mean = ans))
}

