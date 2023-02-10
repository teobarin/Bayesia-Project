library(dplyr)
library(ggplot2)

########## FUNCTIONS #####
f_psi = function(n,k,N_k,m1,lambda, psi, m1_bar){
  S = 0
  for(j in 1:(k-m1)){
    temp = lgamma(N_k[j+m1] - (exp(psi)/(1 + exp(psi)))) - lgamma(1/(1+exp(psi)))
    S = S + temp
  }
  val = dbeta((exp(psi)/(1+exp(psi))),sigma_par$a,sigma_par$b, log=TRUE) + 
    ((k - m1_bar)*(psi - log(1 + exp(psi)))) + S + 
    lgamma((exp(lambda)/exp(psi))+exp(lambda) + k - m1_bar) -
    lgamma((exp(lambda)/exp(psi))+exp(lambda))
  return(val)  
}

g_lambda = function(n,k,N_k,m1,lambda, psi, m1_bar){
  val = dgamma(exp(lambda),theta_par$a,rate=theta_par$b,log=TRUE) + 
    lgamma(exp(lambda)) +
    lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) -
    lgamma((exp(lambda)/exp(psi)) + exp(lambda)) -
    lgamma(exp(lambda) + n - m1_bar)
  return(val)
}

h_m1_bar = function(n,k,N_k,lambda, psi, m1, m1_bar, beta){
  
  v = lchoose(m1, m1_bar) + (n - m1_bar)*log(beta) + m1_bar*log(1-beta) +
    (k-m1_bar)*log(exp(psi)/(1+exp(psi))) +
    lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) - 
    lgamma(exp(lambda) + n - m1_bar)
  
  return(v)
  
}

update_psi = function(n,k,N_k,m1,lambda, psi, m1_bar){
  psi_s = rnorm(1, psi, alpha$psi)
  A = min(1, exp(f_psi(n,k,N_k,m1,lambda, psi_s, m1_bar)  
                 - f_psi(n,k,N_k,m1,lambda, psi, m1_bar) ))
  if(runif(1) < A){
    psi = psi_s
  }
  return(psi)
}

update_lambda = function(n,k,N_k,m1,lambda, psi, m1_bar){
  lambda_s = rnorm(1, lambda, alpha$lambda)
  A = min(1, exp(g_lambda(n,k,N_k,m1,lambda_s, psi, m1_bar) 
                 - g_lambda(n,k,N_k,m1,lambda, psi, m1_bar) ))
  if(runif(1) < A){
    lambda = lambda_s
  }
  return(lambda)
}

update_m1_bar = function(n,k,N_k,lambda, psi, m1, m1_bar, beta){
  X = numeric(m1)
  X = 0*X
  for(i in 1:m1){
    X[i] = h_m1_bar(n,k,N_k,lambda, psi, m1, i, beta)
  }
  s = sample(x = 0:(length(X)-1), size = 1, prob = exp(X - max(X))/sum(exp(X-max(X))))
  
  return (s)
}

f_psi_nc = function(n,k,N_k,m1,lambda,psi){
  S1 = 0
  S2 = 0
  
  for(i in 1:(k-1)){
    temp = log(exp(lambda)+ i*(exp(psi)/(1+exp(psi))))
    S1 = S1 + temp
  }
  
  for(j in 1:(k-m1)){
    temp = lgamma(N_k[j+m1] - (exp(psi)/(1 + exp(psi)))) - lgamma(1/(1+exp(psi)))
    S2 = S2 + temp
  }
  
  val = dbeta((exp(psi)/(1+exp(psi))),sigma_par$a,sigma_par$b, log=TRUE) + S1 + S2
  
  
  return(val)
}

g_lambda_nc = function(n,k,N_k,m1,lambda,psi){
  S=0
  
  for(j in 1:(k-1)){
    temp = log(exp(lambda) + j*(exp(psi)/(1+exp(psi))))
    S = S + temp
  }
  
  val = dgamma(exp(lambda),theta_par$a,rate=theta_par$b,log=TRUE) + 
    lgamma(exp(lambda) + 1) - lgamma(exp(lambda) + n) + S
  
  
  return(val)
}

update_psi_nc = function(n,k,N_k,m1,lambda, psi){
  psi_s = rnorm(1, psi, alpha$psi)
  A = min(1, exp(f_psi_nc(n,k,N_k,m1,lambda, psi_s)  
                 - f_psi_nc(n,k,N_k,m1,lambda, psi) ))
  if(runif(1) < A){
    psi = psi_s
  }
  return(psi)
}

update_lambda_nc= function(n,k,N_k,m1,lambda, psi){
  lambda_s = rnorm(1, lambda, alpha$lambda)
  A = min(1, exp(g_lambda_nc(n,k,N_k,m1,lambda_s, psi) 
                 - g_lambda_nc(n,k,N_k,m1,lambda, psi) ))
  if(runif(1) < A){
    lambda = lambda_s
  }
  return(lambda)
}

cPY = function(data,sigma_par,theta_par,beta_par,alpha,iters){
  
  n = dim(data)[1]
  subset = count(data, data[2])[2]
  ind = which(subset$n==1)
  m1 = length(ind)
  N_k = subset$n
  N_k = sort(N_k)
  k = length(N_k)
  
  sigma = sigma_par$a/(sigma_par$a+sigma_par$b)
  theta = theta_par$a/theta_par$b
  beta = beta_par$a/(beta_par$a+beta_par$b)
  m1_bar = m1*beta
  
  psi = log(sigma) - log(1-sigma)
  lambda = log(theta)
  
  sigma_vec = numeric(iters)
  theta_vec = numeric(iters)
  beta_vec = numeric(iters)
  m1_vec = numeric(iters)
  psi_vec = numeric(iters)
  lambda_vec = numeric(iters)
  
  psi_old = psi
  lambda_old = lambda
  beta_old = beta
  m1_bar_old = m1_bar
  
  ################################################################################
  
  for(i in 1:iters){
    
    psi_new = update_psi(n,k,N_k,m1,lambda_old, psi_old, m1_bar_old)
    sigma_new = exp(psi_new) / (1 + exp(psi_new))
    
    lambda_new = update_lambda(n,k,N_k,m1,lambda_old, psi_new, m1_bar_old)
    theta_new = exp(lambda_new)
    
    beta_new = rbeta(1, beta_par$a + n - m1_bar_old, beta_par$b + m1_bar_old)
    
    m1_bar_new = update_m1_bar(n,k,N_k,lambda_new, psi_new, m1, m1_bar_old, beta_new)
    m1_bar_new = round(m1_bar_new, digits = 0)
    
    psi_vec[i] = psi_new
    lambda_vec[i] = lambda_new
    
    sigma_vec[i] = sigma_new
    theta_vec[i] = theta_new
    beta_vec[i] = beta_new
    m1_vec[i] = m1_bar_new
    
    psi_old = psi_new
    lambda_old = lambda_new
    beta_old = beta_new
    m1_bar_old = m1_bar_new
    
  }
  
  v = list(n=n,N_k=N_k,k=k,m1=m1,sigma=sigma_vec,theta=theta_vec,beta=beta_vec,
           m1_bar=m1_vec,psi=psi_vec,lambda=lambda_vec)
  
  return(v)
  
}

PY = function(data,sigma_par,theta_par,alpha,iters){
  
  n = dim(data)[1]
  subset = count(data, data[2])[2]
  ind = which(subset$n==1)
  m1 = length(ind)
  N_k = subset$n
  N_k = sort(N_k)
  k = length(N_k)
  
  sigma = sigma_par$a/(sigma_par$a+sigma_par$b)
  theta = theta_par$a/theta_par$b
  
  psi = log(sigma) - log(1-sigma)
  lambda = log(theta)
  
  sigma_vec_nc = numeric(iters)
  theta_vec_nc = numeric(iters)
  psi_vec_nc = numeric(iters)
  lambda_vec_nc = numeric(iters)
  
  psi_old_nc = psi
  lambda_old_nc = lambda
  
  ################################################################################
  
  for(i in 1:iters){
    
    psi_new_nc = update_psi_nc(n,k,N_k,m1,lambda_old_nc, psi_old_nc)
    sigma_new_nc = exp(psi_new_nc) / (1 + exp(psi_new_nc))
    
    lambda_new_nc = update_lambda_nc(n,k,N_k,m1,lambda_old_nc, psi_new_nc)
    theta_new_nc = exp(lambda_new_nc)
    
    sigma_vec_nc[i] = sigma_new_nc
    theta_vec_nc[i] = theta_new_nc
    
    psi_vec_nc[i] = psi_new_nc
    lambda_vec_nc[i] = lambda_new_nc
    
    psi_old_nc = psi_new_nc
    lambda_old_nc = lambda_new_nc
    
  }
  
  v = list(n=n,N_k=N_k,k=k,m1=m1,sigma_nc=sigma_vec_nc,theta_nc=theta_vec_nc,
           psi_nc=psi_vec_nc,lambda_nc=lambda_vec_nc)
  
  return(v)
  
}

expected_m_py <- Vectorize(function(m, n, sigma, alpha) {
  out <- log(alpha) + lchoose(n, m) + lgamma(m - sigma) - lgamma(1 - sigma) - lgamma(alpha + n) + lgamma(alpha) + lgamma(alpha + sigma + n - m) - lgamma(alpha + sigma)
  exp(out)
}, vectorize.args = "m")

# Importa dataset
data <- read.csv("india_small.csv", sep ="\t")
#names <- (data$verbatimScientificName)
first_column <- c('')
second_column = data['verbatimScientificName'] 
df <- data.frame(first_column, second_column)

sigma_par = list(a=1,b=1)
theta_par = list(a=2,b=0.02)
beta_par = list(a=1,b=1)
alpha = list(psi = 0.2, lambda = 0.4)
iters = 5000

d = cPY(df,sigma_par,theta_par,beta_par,alpha,iters)
d_nc = PY(df,sigma_par,theta_par,alpha,iters)


#grafico log-log
freq = d$N_k
n = sum(freq)
M_l <- as.numeric(table(factor(freq, levels = 1:n))) #numero di cluster con frequenza l=1,...,n
idx <- 1:(which.min(M_l) + 30) # which(M_l > 0)
col <- c("#009E73", "#D55E00", "#000000", "#0072B2")

est_df_temp <- data.frame(
  x = idx, 
  y = expected_m_py(idx, n = sum(freq), sigma = quantile(d$sigma, 0.5), alpha = quantile(d$theta, 0.5)), 
  y_low = expected_m_py(idx, n = sum(freq), sigma = quantile(d$sigma, 0.05), alpha = quantile(d$theta, 0.05)), 
  y_up = expected_m_py(idx, n = sum(freq), sigma = quantile(d$sigma, 0.95), alpha = quantile(d$theta, 0.95)),
  yPY = expected_m_py(idx, n = sum(freq), sigma = quantile(d_nc$sigma_nc, 0.5), alpha = quantile(d_nc$theta_nc, 0.5)), 
  y_lowPY = expected_m_py(idx, n = sum(freq), sigma = quantile(d_nc$sigma_nc, 0.05), alpha = quantile(d_nc$theta_nc, 0.05)), 
  y_upPY = expected_m_py(idx, n = sum(freq), sigma = quantile(d_nc$sigma_nc, 0.95), alpha = quantile(d_nc$theta_nc, 0.95))
)

ggplot(data.frame(x = idx[M_l != 0], y = M_l[idx[M_l != 0]]), aes(x = x, y = y)) +  
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_line(data = est_df_temp,
            mapping = aes(x = x, y = yPY), col = col[2]) +
  geom_ribbon(data = est_df_temp,
              mapping = aes(x = x, ymin = y_lowPY, ymax = y_upPY), fill = col[2], alpha = 0.2) +
  geom_line(data = est_df_temp,
            mapping = aes(x = x, y = y), col = col[1]) +
  geom_ribbon(data = est_df_temp,
              mapping = aes(x = x, ymin = y_low, ymax = y_up), fill = col[1], alpha = 0.2) +
  geom_segment(data = est_df_temp, mapping = aes(x = 1, xend = 2, y = y[1] + mean(d$m1_bar) , 
                                                 yend = y[2]), col = col[1], lty = 2) +
  geom_point(shape = 1, alpha = 0.9, col = col[3]) + 
  geom_point(shape = 16, alpha = 0.3, col = col[3]) +
  xlab(expression(n[j] * ' (log-scale)')) + 
  ylab(expression(m[j] * ' (log-scale)'))
