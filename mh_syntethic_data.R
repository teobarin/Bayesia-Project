library(dplyr)

## Analizzo il dataset e recupero i valori n, k, m_1
data = read.csv(paste0("synthetic_new_", 1,".csv"), header = TRUE)
n = dim(data)[1]

subset = count(data,data[2])[2]
ind = which(subset$n==1)
m1 = length(ind)

N_k = subset$n
N_k = sort(N_k)
k = length(N_k)

################################################################################
## ALGORITHM 1
## Definisco i parametri delle prior
## psi = log(sigma/(1-sigma))
## lambda = log(theta)

theta <- 100
sigma <- 0.2
beta <- 0.95

a <- list(theta = 2, sigma = 1, beta = 1)
b <- list(theta = 0.02, sigma = 1, beta = 1)

psi <- log(sigma) - log(1 - sigma)
lambda <- log(theta)

alpha <- list(theta = 0.4, sigma = 0.2)
#alpha <- list(theta = 1, sigma = 0.5)

#m1_bar = round(runif(1,0,m1), digits=0)
m1_bar = round((1-beta)*n)


#### FUNCTIONS

f_psi = function(lambda, psi, m1_bar){
  
  S = 0
  for(j in 1:(k-m1)){
    temp = lgamma(N_k[j+m1] - (exp(psi)/(1 + exp(psi)))) - lgamma(1/(1+exp(psi)))
    S = S + temp
  }
  
  val = dbeta((exp(psi)/(1+exp(psi))),a$sigma,b$sigma, log=TRUE) + 
        ((k - m1_bar)*(psi - log(1 + exp(psi)))) + S + 
        lgamma((exp(lambda)/exp(psi))+exp(lambda) + k - m1_bar) -
        lgamma((exp(lambda)/exp(psi))+exp(lambda))
  
  return(val)  
}

g_lambda = function(lambda, psi, m1_bar){
  
  val = dgamma(exp(lambda),a$theta,rate=b$theta,log=TRUE) + 
        lgamma(exp(lambda)) +
        lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) -
        lgamma((exp(lambda)/exp(psi)) + exp(lambda)) -
        lgamma(exp(lambda) + n - m1_bar)
  
  return(val)
}

update_psi = function(lambda, psi, m1_bar, accept){
  psi_s = rnorm(1, psi, alpha$sigma)
  
  A = min(1,exp(f_psi(lambda,psi_s,m1_bar)  
                - f_psi(lambda,psi,m1_bar) ))
  if(runif(1) < A){
    psi = psi_s
    accept = accept + 1
  }
  return(list(a = psi, b = accept))
}

update_lambda = function(lambda, psi, m1_bar, accept){

  lambda_s = rnorm(1, lambda, alpha$theta)
  
  A = min(1, exp(g_lambda(lambda_s, psi, m1_bar) 
                - g_lambda(lambda, psi, m1_bar) ))
  if(runif(1) < A){
    lambda = lambda_s
    accept = accept + 1
  }
  return(list(a = lambda, b = accept))
}

######## FUNZIONE M1_BAR ########
h  = function(lambda, psi, m1, m1_bar, beta){
  
  v = lchoose(m1, m1_bar) + (n - m1_bar)*log(beta) + m1_bar*log(1-beta) +
    (k-m1_bar)*log(exp(psi)/(1+exp(psi))) +
    lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) - 
    lgamma(exp(lambda) + n - m1_bar)
  
  return(v)
  
}

up_m1_bar = function(lambda, psi, m1, m1_bar, beta){
  
  X = 0*X
  for(i in 1:m1){
    X[i] = h(lambda, psi, m1, i, beta)
  }
  s = sample(x = 0:(length(X)-1), size = 1, prob = exp(X - max(X))/sum(exp(X-max(X))))
  
  return (s)
}


iters = 8000

sigma_vec = numeric(iters)
theta_vec = numeric(iters)
beta_vec = numeric(iters)
m1_vec = numeric(iters)
psi_vec = numeric(iters)
lambda_vec = numeric(iters)
X = numeric(m1)

psi_old = psi
lambda_old = lambda
beta_old = beta
m1_bar_old = m1_bar

accept_theta = 0
accept_sigma = 0

for(i in 1:iters){
  
  new = list(a, b)
  new = update_psi(lambda_old, psi_old, m1_bar_old, accept_theta)
  psi_new = new$a
  accept_theta = new$b
  sigma_new = exp(psi_new)/(1+exp(psi_new))
  
  new = list(a, b)
  new = update_lambda(lambda_old, psi_new, m1_bar_old, accept_sigma)
  lambda_new = new$a
  accept_sigma = new$b
  theta_new = exp(lambda_new)
  
  #beta_new = exp(dbeta(beta_old,a$beta + n - m1_bar, b$beta + m1_bar_old, log=TRUE) + lbeta(a$beta + n - m1_bar_old,b$beta + m1_bar_old)-lbeta(a$beta,b$beta))
  #beta_new = dbeta(beta_old, a$beta + n - m1_bar_old, b$beta + m1_bar_old)
  #beta_new = beta_old
  beta_new = rbeta(1, a$beta + n - m1_bar_old, b$beta + m1_bar_old)
  
  #m1_bar_new = m1_bar_old
  m1_bar_new = up_m1_bar(lambda_new, psi_new, m1, m1_bar_old, beta_new)
  m1_bar_new = round(m1_bar_new, digits = 0)
  
  sigma_vec[i] = sigma_new
  theta_vec[i] = theta_new
  psi_vec[i] = psi_new
  lambda_vec[i] = lambda_new
  beta_vec[i] = beta_new
  m1_vec[i] = m1_bar_new
  
  psi_old = psi_new
  lambda_old = lambda_new
  beta_old = beta_new
  m1_bar_old = m1_bar_new
  
}

print(accept_theta/iters)
print(accept_sigma/iters)


par(mfrow=c(2,2))
plot(sigma_vec, main=paste("sigma value, mean = ", mean(sigma_vec)), type="l")
plot(theta_vec, main=paste("theta value, mean = ", mean(theta_vec)), type="l" )
plot(psi_vec, main=paste("psi value, mean = ", mean(psi_vec)), type="l")
plot(lambda_vec, main=paste("lambda value, mean = ", mean(lambda_vec)), type="l" )
#plot(beta_vec, main=paste("beta value, mean = ", mean(beta_vec)), type="l")
#plot(m1_vec, main=paste("m1_bar value, mean = ", mean(m1_vec)), type="l")

mean(sigma_vec)
mean(theta_vec)
mean(beta_vec)
mean(m1_vec)