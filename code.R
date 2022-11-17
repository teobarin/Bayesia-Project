library(readxl)
library(catdata)
library(dplyr)
data = read_excel("Dataset.xlsx")
n=dim(data)[1]

subset=count(data,verbatimScientificName)
ind=which(subset$n==1)
m1=length(ind)

N_k=subset[[2]]
N_k=sort(N_k)
k=length(N_k)

################################################################################
# Algoritmo su R cos? da non dover modificare MH in Python

## ALGORITHM 1
## MH PER UPDATE DI THETA E SIGMA
## psi = log(sigma/(1-sigma))
## lambda = log(theta)

a_theta <- 2
b_theta <- 0.02
a_sigma <- 1
b_sigma <- 1
theta <- 100
lambda <- log(theta)
sigma <- 0.6
n <- 50000
k <- 619
m1 <- 500
beta <- 0.5
m1_bar <- runif(1,0,m1)
a_beta <- 1
b_beta <- 1
#a_beta_post <- a_beta_prior+n-m1_bar
#b_beta_post <- b_beta_prior+m1_bar
alpha <- 10


## UPDATE SIGMA
f_sigma = function(sigma, theta, m1_bar){
  psi <- log(sigma/(1-sigma))
  S = 0
  for (i in (m1+1):k){
    S = S + lgamma(N_k[i] - sigma)
  }
  res = (a_sigma+k-m1_bar-1)*psi - (a_sigma+b_sigma+k-m1_bar-2)*log(1+exp(psi)) +
      lgamma(a_sigma+b_sigma) - lgamma(a_sigma) - lgamma(b_sigma) + S -
      (k-m1-1)*lgamma(1-sigma) + lgamma(theta/sigma+k-m1_bar) - lgamma(theta/sigma)
  res
}

mh_sigma=function(iters, alpha, sigma, theta, m1_bar)
{
  xvec=numeric(iters)

  for (i in 1:iters) {
    sigma_s = sigma + runif(1,0,alpha)
    A = f_sigma(sigma_s, theta, m1_bar)/f_sigma(sigma, theta, m1_bar)
    if (runif(1) < A)
      sigma = sigma_s
    xvec[i] = sigma
  }
  
  return(xvec)
}

post_s = mh_sigma(n, alpha, sigma, theta, m1_bar)
hist(post_s,100,freq=FALSE, main="mh_sigma")
curve(dbeta(x,a_sigma,b_sigma),add=TRUE,col=2,lwd=2)


## UPDATE THETA

f_theta = function(x, lambda, psi, m1_bar){
  num = lgamma(lambda)+lgamma((lambda/psi)+k-m1_bar)
  den = lgamma(lambda/psi)+lgamma(lambda+n-m1_bar)
  res = dgamma(x,a_theta,b_theta)*(num-den)
  res
}

mh_theta = function(iters,alpha, lambda, psi, m1_bar)
{
  xvec = numeric(iters)
  x = 1
  for (i in 1:iters) {
    xs = x+rnorm(1,0,alpha)
    A = f_theta(xs, lambda, psi, m1_bar)/f_theta(x, lambda, psi, m1_bar)
    if (runif(1)<A)
      x = xs
    xvec[i] = x
  }
  return(xvec)
}

#post_t = mh_theta(n,alpha, lambda, psi, m1_bar)
#hist(post_t,100,freq=FALSE,main="mh_theta")
#curve(dgamma(x,a_theta,b_theta),add=TRUE,col=2,lwd=2)


# beta e m1 sampliamo solo un valore random?
R = 1

### valori iniziali
lambda_old = lambda
psi_old = psi
beta_old = beta
m1_bar_old = m1_bar

for (r in 1:R){
  
  # update sigma
  psi_new = mh_sigma(n, alpha, lambda_old, psi_old, m1_bar_old)
  psi_new = mean(psi_new)
  
  # NB: rising factorial -> (1-sigma)_(1-n(i)) = gamma(1-n(i))/(1-sigma)
  # chi sono questi benedetti n(i)??
  
  # update theta
  lambda_new = mh_theta(n, alpha, lambda_old, psi_new, m1_bar_old)
  lambda_new = mean(lambda_new)
  
  # update beta -> Expected values of beta
  #beta_new = (a_beta + n - m1_bar_old)/(a_beta + n + b_beta)
  beta_new = beta_old
  
  # update m1_bar
  num = lgamma((lambda_new/psi_new) + k - m1_bar_old)
  den = lgamma(lambda_new + n - m1_bar_old)
  
  m1_bar_old = round(m1_bar_old, digits = 0)
  
  m1_bar_new = log(choose(m1,m1_bar_old)) + log(beta_new)*(n-m1_bar_old) + log((1-beta_new))*(m1_bar_old) + psi_new*(k-m1_bar_old) + num - den
  m1_bar_new = exp(m1_bar_new)
  #m1_bar_new = (choose(m1,m1_bar_old))*(beta_new^(n-m1_bar_old))*((1-beta_new)^(m1_bar_old))*(exp(psi_new)^(k-m1_bar_old))*num/den
  
  lambda_old <- lambda_new
  psi_old <- psi_new
  beta_old <- beta_new
  m1_bar_old <- m1_bar_new
}


