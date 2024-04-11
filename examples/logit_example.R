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

# Ground 'Truth'
print(colMeans(chain_steps))
# 'Random Samples from Chain Output'
print(colMeans(theta))
# Apply importance sampling
out = BBIS(theta,theta_grads,1000)
print(out$adj_mean)
