library(cmdstanr)
library(mcmc)
library(GFBBIS)


data(logit)

set.seed(12345)

chain_len = 100000

mod<-cmdstan_model("~/GFBBIS/examples/logit_example.stan", compile_model_methods=TRUE,force_recompile = TRUE)

chain<-mod$sample(data=list(Y=logit[,1],X=logit[,2:5]), chains=1, seed=12345, iter_sampling=chain_len, iter_warmup=10000)

chain$init_model_methods()

chain_steps = matrix(data=NA,nrow=chain_len,ncol=5)
chain_steps[,1]=chain$draws('y_intercept')
chain_steps[,2]=chain$draws('beta[1]')
chain_steps[,3]=chain$draws('beta[2]')
chain_steps[,4]=chain$draws('beta[3]')
chain_steps[,5]=chain$draws('beta[4]')

grad_func <- function(x0) {
  output_matrix = matrix(data=NA, nrow=nrow(x0),ncol=ncol(x0))

  for (ii in 1:nrow(x0)) {
    output_matrix[ii,] = chain$grad_log_prob(x0[ii,])
  }

  return(output_matrix)
}

chain_grads <- grad_func(chain_steps)

iis = sample(1:nrow(chain_steps),50)

theta = chain_steps[iis,]
theta_grads = chain_grads[iis,]

# Base BBIS Example
# Ground 'Truth'
print(colMeans(chain_steps))
# 'Random Samples from Chain Output'
print(colMeans(theta))
# Apply importance sampling
out = BBIS(theta,theta_grads,1000)
print(out$adj_mean)


mvtnorm_prob <- function(theta,mu,sigma) {
  probs = matrix(data=NA,nrow=nrow(theta),ncol=1)
  grad_probs = matrix(data=NA,nrow=nrow(theta),ncol=ncol(theta))


  for (ii in 1:nrow(theta)) {
    probs[ii] = ((2*pi)**(-ncol(theta)/2)) * (det(sigma)**(-1/2)) * exp((-1/2)*(t(theta[ii,]-mu)%*%solve(sigma)%*%(theta[ii,]-mu)))

    grad_probs[ii,] = -solve(sigma)%*%(theta[ii,]-mu)
  }

  return(list(den=(probs),grad_x=grad_probs))
}

# Gradient Free BBIS Example, using surrogate
iis2 = sample(1:nrow(chain_steps),50)
theta2 = chain_steps[iis2,]
theta_probs= mvtnorm_prob(theta2, mu=colMeans(chain_steps),sigma=var(chain_steps))

out2 = GFBBIS(theta2,theta_probs$den,1000)
print(colMeans(chain_steps))
print(colMeans(theta2))
print(out2$adj_mean)





# Gradient Free BBIS Example, no surrogate
iis3 = sample(1:nrow(chain_steps),50)
theta3 = chain_steps[iis3,]*1
theta_probs = mvtnorm_prob(theta3,mu=colMeans(chain_steps),sigma=var(chain_steps))

out3 = FGFBBIS(theta3,theta_probs$den,1000)
print(colMeans(chain_steps))
print(colMeans(theta3))
print(out3$adj_mean)













MSEs_standard = matrix(data=NA,nrow=500,ncol=1)
MSEs_base = matrix(data=NA,nrow=500,ncol=1)
MSEs_gf = matrix(data=NA,nrow=500,ncol=1)
MSEs_fgf = matrix(data=NA,nrow=500,ncol=1)

for (ii in 1:500) {
  print(ii)

  iis = sample(1:nrow(chain_steps),50)

  theta = chain_steps[iis,]*1.5
  theta_grads = grad_func(theta)
  theta_probs = mvtnorm_prob(theta,mu=colMeans(chain_steps),sigma=var(chain_steps))

  out_base = BBIS(theta,theta_grads,1000)
  out_gf = GFBBIS(theta,theta_probs$den,1000)
  out_fgf = FGFBBIS(theta,theta_probs$den,1000)

  MSEs_standard[ii,] = sum((colMeans(theta)-colMeans(chain_steps))**2)
  MSEs_base[ii,] = sum((out_base$adj_mean-colMeans(chain_steps))**2)
  MSEs_gf[ii,] = sum((out_gf$adj_mean-colMeans(chain_steps))**2)
  MSEs_fgf[ii,] = sum((out_fgf$adj_mean-colMeans(chain_steps))**2)
}


plot(1:500,log(MSEs_standard),ylim = c(-7,0),col='black',ylab = 'log MSEs')
abline(h=mean(log(MSEs_standard)),col='black')
points(log(MSEs_base),col='blue')
abline(h=mean(log(MSEs_base)),col='blue')
points(log(MSEs_gf),col='red')
abline(h=mean(log(MSEs_gf)),col='red')
points(log(MSEs_fgf),col='purple')
abline(h=mean(log(MSEs_fgf)),col='purple')
legend(1,-5,legend = c('No IS', 'Base IS', 'GF IS', 'FGF IS'), col=c('black','blue','red','purple'),lty=1)






