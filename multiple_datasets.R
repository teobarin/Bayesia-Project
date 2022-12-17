library(readxl)
library(catdata)
library(dplyr)
library(ggplot2)
library(corrplot)

## ALGORITHM 1
## Definisco i parametri delle prior
## psi = log(sigma/(1-sigma))
## lambda = log(theta)

theta <- 100
sigma <- 0.2
beta <- 0.95

a <- list(theta = 2, sigma = 1, beta = 1)
b <- list(theta = 0.02, sigma = 5, beta = 1)

psi <- log(sigma) - log(1-sigma)
lambda <- log(theta)
alpha <- list(theta = 1, sigma = 0.5)

f_psi = function(lambda,psi,m1_bar){
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

g_lambda = function(lambda,psi,m1_bar){
  val = dgamma(exp(lambda),a$theta,rate=b$theta,log=TRUE) + 
    lgamma(exp(lambda)) +
    lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) -
    lgamma((exp(lambda)/exp(psi)) + exp(lambda)) -
    lgamma(exp(lambda) + n - m1_bar)
  return(val)
}

update_psi = function(lambda,psi,m1_bar){
  psi_s = rnorm(1,psi,alpha$sigma)
  A = min(1,exp(f_psi(lambda,psi_s,m1_bar)  
                - f_psi(lambda,psi,m1_bar) ))
  if(runif(1) < A){
    psi = psi_s
  }
  return(psi)
}

update_lambda = function(lambda,psi,m1_bar){
  lambda_s = rnorm(1,lambda,alpha$theta)
  A = min(1,exp(g_lambda(lambda_s,psi,m1_bar) 
                - g_lambda(lambda,psi,m1_bar) ))
  if(runif(1) < A){
    lambda = lambda_s
  }
  return(lambda)
}

sigma_values = numeric(10)
theta_values = numeric(10)

sigma_list = list()
theta_list = list()

## Analizzo il dataset e recupero i valori n, k, m_1

for(j in 0:9){
  
  data = read.csv(paste0("synthetic_data_", j,".csv"), header = TRUE)
  n = dim(data)[1]
  
  subset = count(data,data[2])[2]
  ind = which(subset$n==1)
  m1 = length(ind)
  m1_bar = round(runif(1,0,m1), digits=0)
  
  N_k = subset$n
  N_k = sort(N_k)
  k = length(N_k)
  
  iters = 5000
  
  sigma_vec = numeric(iters)
  theta_vec = numeric(iters)
  beta_vec = numeric(iters)
  m1_vec = numeric(iters)
  
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
    #beta_new = beta_old
    
    #beta_new = rbeta(1, a$beta + n - m1_bar_old, b$beta + m1_bar_old)
    
    #m1_bar_new = exp(lchoose(m1,m1_bar_old) + (log(beta_new)*(n-m1_bar_old)) + (log(1 - beta_new)*(m1_bar_old)) + (log(sigma_new)*(k-m1_bar_old)) + lgamma((theta_new/sigma_new) + k - m1_bar_old) - lgamma(theta_new + n - m1_bar_old))
    #m1_bar_new = m1_bar_old
    
    psi_old = psi_new
    lambda_old = lambda_new
    
    #beta_old = beta_new
    #m1_bar_old = round(m1_bar_new, digits=0)
    
    sigma_vec[i] = sigma_new
    theta_vec[i] = theta_new
    #beta_vec[i] = beta_new
    #m1_vec[i] = m1_bar_new
  }
  
  sigma_list[j+1] = list(c(sigma_vec))
  theta_list[j+1] = list(c(theta_vec))

  sigma_values[j+1] = mean(sigma_vec)
  theta_values[j+1] = mean(theta_vec)
  
}

sigma_values
theta_values

################################################################################

par(mfrow=c(2,5))

for(j in 0:9){
  
  plot(sigma_list[[j+1]], main="sigma value", type="l")
  #plot(theta_list[[j+1]], main= "theta value", type="l")
  
}


plot(sigma_list[[4+1]], main="sigma value", type="l")
