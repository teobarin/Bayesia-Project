# psi = log(sigma)-log(1-sigma)
# lambda = log(theta)


### updating function for psi in contaminated model
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


### updating function for lambda in contaminated model
g_lambda = function(n,k,N_k,m1,lambda, psi, m1_bar){
  val = dgamma(exp(lambda),theta_par$a,rate=theta_par$b,log=TRUE) + 
    lgamma(exp(lambda)) +
    lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) -
    lgamma((exp(lambda)/exp(psi)) + exp(lambda)) -
    lgamma(exp(lambda) + n - m1_bar)
  return(val)
}


### updating function for m1_bar
h_m1_bar = function(n,k,N_k,lambda, psi, m1, m1_bar, beta){  
  v = lchoose(m1, m1_bar) + (n - m1_bar)*log(beta) + m1_bar*log(1-beta) +
    (k-m1_bar)*log(exp(psi)/(1+exp(psi))) +
    lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) - 
    lgamma(exp(lambda) + n - m1_bar)
  return(v)
}


### Metropolis-Hastings for psi in contaminated model
update_psi = function(n,k,N_k,m1,lambda, psi, m1_bar){
  psi_s = rnorm(1, psi, alpha$psi)
  A = min(1, exp(f_psi(n,k,N_k,m1,lambda, psi_s, m1_bar)  
                 - f_psi(n,k,N_k,m1,lambda, psi, m1_bar) ))
  if(runif(1) < A){
    psi = psi_s
  }
  return(psi)
}


### Metropolis-Hastings for lambda in contaminated model
update_lambda = function(n,k,N_k,m1,lambda, psi, m1_bar){
  lambda_s = rnorm(1, lambda, alpha$lambda)
  A = min(1, exp(g_lambda(n,k,N_k,m1,lambda_s, psi, m1_bar) 
                 - g_lambda(n,k,N_k,m1,lambda, psi, m1_bar) ))
  if(runif(1) < A){
    lambda = lambda_s
  }
  return(lambda)
}


### Function that sample m1_bar
update_m1_bar = function(n,k,N_k,lambda, psi, m1, m1_bar, beta){
  X = numeric(m1)
  X = 0*X
  for(i in 1:m1){
    X[i] = h_m1_bar(n,k,N_k,lambda, psi, m1, i, beta)
  }
  s = sample(x = 0:(length(X)-1), size = 1, prob = exp(X - max(X))/sum(exp(X-max(X))))
  return (s)
}


### updating function for psi in non contaminated model
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


### updating function for lambda in non contaminated model
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


### Metropolis-Hastings for psi in non contaminated model
update_psi_nc = function(n,k,N_k,m1,lambda, psi){
  psi_s = rnorm(1, psi, alpha$psi)
  A = min(1, exp(f_psi_nc(n,k,N_k,m1,lambda, psi_s)  
                 - f_psi_nc(n,k,N_k,m1,lambda, psi) ))
  if(runif(1) < A){
    psi = psi_s
  }
  return(psi)
}


### Metropolis-Hastings for psi in non contaminated model
update_lambda_nc= function(n,k,N_k,m1,lambda, psi){
  lambda_s = rnorm(1, lambda, alpha$lambda)
  A = min(1, exp(g_lambda_nc(n,k,N_k,m1,lambda_s, psi) 
                 - g_lambda_nc(n,k,N_k,m1,lambda, psi) ))
  if(runif(1) < A){
    lambda = lambda_s
  }
  return(lambda)
}


### Function that from the frequencies of the clusters recreate a vector that represent
### the data and return the data
extract_data = function(freq){  
  l = length(freq)
  x = seq(l)
  x = as.factor(x)
  data = as.data.frame(cbind(x,freq))
  new_data = seq(l)
  for(i in 1:l){
    f = data$freq[i]
    v = numeric(f-1) + i
    new_data = c(new_data,v)
  }
  new_data = sort(new_data)
  new_data = as.factor(new_data)
  new_data = as.data.frame(new_data)
  return(new_data)
}


### Create "it" different frequencies from a dataset and return the list
### that contain all the frequencies
create_freq = function(it,data,m){
  empty_list <- vector("list", length = 5)
  n = dim(data)[1]
  for(i in 1:it){
    ind = sample(1:n,m)
    g = as.data.frame(data[-ind,])
    names(g)=c("dati")
    sub = count(g,dati)[2]
    N_k = sub$n
    empty_list[[i]] <- N_k
  }
  return(empty_list)
}


### Function that given parameters estimate with the model returned
### the prediction of the number of singleton in the next m data
pred_singl = function(n,m,sigma,theta,beta,k,M_n){
  S = 0
  p = log(theta+(k-M_n)*sigma) + lgamma(theta+n-M_n) - lgamma(theta+n-M_n+sigma)
  for (l in 0:m){
    a = lchoose(m,l) + l*log(1-beta) + (m-l)*log(beta)
    b = log(m-l) +  lgamma(theta+n-M_n+sigma+m-l-1)   - lgamma(theta+n-M_n+m-l) + p
    t = (exp(b) + l) * exp(a)
    S = S + t
  }
  return (round(S,digits=0))
}


### Function that given parameters estimate with the model returned
### the prediction of the number of clusters with r elements in the next m data
pred_dim_clust = function(n,m,r,sigma,theta,beta,k,M_n){
  S=0
  p = lgamma(r-sigma) - lgamma(1-sigma) +log(theta+(k-M_n)*sigma) 
  + lgamma(theta+n-M_n) - lgamma(theta+n-M_n+sigma)
  for(l in 0:(m-r)){
    a = lchoose(m,l) + l*log(1-beta) + (m-l)*log(beta)
    b = lchoose(m-l,r)  +  lgamma(theta+n-M_n+sigma+m-l-r) - lgamma(theta+n-M_n+m-l) + p
    t = exp(a + b)
    S = S + t
  }
  return(round(S,digits=0))
} 


### Function that given parameters estimate with the model returned
### the prediction of the number of clusters in the next m data
pred_clust = function(n,m,sigma,theta,beta,k,M_n){
  S = 0
  c = k - M_n + theta/sigma
  p = log(c) + lgamma(theta+n-M_n) - lgamma(theta+n-M_n+sigma)
  for (l in 0:m){
    a = lchoose(m,l) + l*log(1-beta) + (m-l)*log(beta)
    b = lgamma(theta+n-M_n+sigma+m-l) - lgamma(theta+n-M_n+m-l) + p
    t = (exp(b) + l - c) * exp(a)
    S = S + t
  }
  return (round(S,digits=0))
}


### Function that give the frequencies of the clusters returns the vectors 
### with all values of the parameters sampled for the contaminated model
cPY_freq = function(freq,sigma_par,theta_par,beta_par,alpha,iters){ 
  n = sum(freq)
  ind = which(freq==1)
  m1 = length(ind)
  N_k = freq
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


### Function that give the frequencies of the clusters returns the vectors 
### with all values of the parameters sampled for the non contaminated model
PY_freq = function(freq,sigma_par,theta_par,alpha,iters){  
  n = sum(freq)
  ind = which(freq==1)
  m1 = length(ind)
  N_k = freq
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

