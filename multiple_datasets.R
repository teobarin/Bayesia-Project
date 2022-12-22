library(dplyr)

## ALGORITHM 1
## Definisco i parametri delle prior
## psi = log(sigma/(1-sigma))
## lambda = log(theta)

theta <- 100
sigma <- 0.2
beta <- 0.95

a <- list(theta = 2, sigma = 1, beta = 1)
b <- list(theta = 0.02, sigma = 1, beta = 1)

psi <- log(sigma) - log(1-sigma)
lambda <- log(theta)
alpha <- list(theta = 0.4, sigma = 0.2)

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

update_psi = function(lambda, psi, m1_bar){
  psi_s = rnorm(1, psi, alpha$sigma)
  A = min(1, exp(f_psi(lambda, psi_s, m1_bar)  
                - f_psi(lambda, psi, m1_bar) ))
  if(runif(1) < A){
    psi = psi_s
  }
  return(psi)
}

update_lambda = function(lambda, psi, m1_bar){
  lambda_s = rnorm(1, lambda, alpha$theta)
  A = min(1, exp(g_lambda(lambda_s, psi, m1_bar) 
                - g_lambda(lambda, psi, m1_bar) ))
  if(runif(1) < A){
    lambda = lambda_s
  }
  return(lambda)
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


sigma_values = numeric(10)
theta_values = numeric(10)
beta_values = numeric(10)
m1_bar_values = numeric(10)

sigma_list = list()
theta_list = list()
beta_list = list()
m1_bar_list = list()

## Analizzo il dataset e recupero i valori n, k, m_1

for(j in 0:9){
  
  data = read.csv(paste0("synthetic_new_", j,".csv"), header = TRUE)
  n = dim(data)[1]
  
  subset = count(data, data[2])[2]
  ind = which(subset$n==1)
  m1 = length(ind)
  m1_bar = round(runif(1, 0, m1), digits=0)
  
  N_k = subset$n
  N_k = sort(N_k)
  k = length(N_k)
  
  iters = 5000
  
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
  
  ################################################################################
  
  for(i in 1:iters){
    
    psi_new = update_psi(lambda_old, psi_old, m1_bar_old)
    sigma_new = exp(psi_new) / (1 + exp(psi_new))
    
    lambda_new = update_lambda(lambda_old, psi_new, m1_bar_old)
    theta_new = exp(lambda_new)
    
    #beta_new = exp(dbeta(beta_old,a$beta + n - m1_bar, b$beta + m1_bar_old, log=TRUE) + lbeta(a$beta + n - m1_bar_old,b$beta + m1_bar_old)-lbeta(a$beta,b$beta))
    #beta_new = dbeta(beta_old,a$beta + n - m1_bar_old, b$beta + m1_bar_old)
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
  
  sigma_list[j+1] = list(c(sigma_vec))
  theta_list[j+1] = list(c(theta_vec))
  beta_list[j+1] = list(c(beta_vec))
  m1_bar_list[j+1] = list(c(m1_vec))

  sigma_values[j+1] = mean(sigma_vec)
  theta_values[j+1] = mean(theta_vec)
  beta_values[j+1] = mean(beta_vec)
  m1_bar_values[j+1] = mean(m1_bar_vec)
  
}

sigma_values
theta_values

################################################################################

par(mfrow=c(2,5))

for(j in 0:9){
  
  #plot(sigma_list[[j+1]], main=paste("sigma", j, "-> mean = ", signif(sigma_values[j+1], digits = 3)), type="l")
  plot(theta_list[[j+1]], main=paste("theta", j, " -> mean = ", trunc(theta_values[j+1])), type="l")
  #plot(beta_list[[j+1]], main=paste("beta", j, " -> mean = ", signif(beta_values[j+1], digits = 3)), type="l")
  #plot(m1_bar_list[[j+1]], main=paste("m1_bar", j, " -> mean = ", trunc(m1_bar_values[j+1])), type="l")
  
}

