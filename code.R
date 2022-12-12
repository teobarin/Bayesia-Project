#
# SIMBOLI: α β λ θ ψ σ
# 
# λ = log(θ)     θ = exp(λ)
#                                       exp(ψ)
# ψ = log(σ) - log(1-σ)         σ = -------------
#                                     1 + exp(ψ)
#
############################# METROPOLIS HASTINGS ##############################
#
# 1) k=1, define x_1
# 2) find x* from q(x*|x_k) and calculate
#
#                         π(x*)q(x_k|x*)    
#     α(x_k,x*)= min ( 1, -------------- )
#                         π(x)q(x*|x_k)
#
# 3) accept x* with probability α(x_k,x*) using u ~ U(0,1):
#    if α(x_k,x*) >= u, x_k+1 = x*
#    if α(x_k,x*) < u, x_k+1 = x_k
#
# 4) set k = k+1 and go back to step 2)
#
################################# UPDATE SIGMA #################################
#
# prior density σ: σ ~ Beta(1,1)
#
# proposal density q(x*|x) ~ N(x,α)
#
#                                  k                         Г(θ/σ + k - m1_bar)
# f(σ) = L(σ) * σ^(k-m1_bar) *     ∏    (1 - σ)_(n_i - 1) * ---------------------
#                             i=m1_bar+1                           Г(θ/σ)
#
# posterior density σ: L(σ|X,θ,m1_bar) ~ L(σ)*f(σ)
#
# posterior density λ: log(L(σ|X,θ,m1_bar)) ~ log(L(exp(ψ)/(1 + exp(ψ))) + log(f(exp(ψ)/(1 + exp(ψ))))
#
################################# UPDATE THETA #################################
#
# prior density θ: θ ~ Gamma(2,0.02)
#
# proposal density q(x*|x) ~ N(x,α)
#
#         Г(θ) * Г(θ/σ + k - m1_bar)
# g(θ) = ----------------------------
#         Г(θ/σ) * Г(θ + n - m1_bar)
#
# posterior density θ: L(θ|X,σ,m1_bar) ~ L(θ)*g(θ)
#
# posterior density λ: log(L(θ|X,σ,m1_bar)) ~ log(L(exp(λ))) + log(g(exp(λ)))
#
################################################################################
# RECUPERO n,k,m1_bar DA DATASET
# data = read_excel("Dataset.xlsx")
# n = dim(data)[1]
# subset = count(data,verbatimScientificName)
# ind = which(subset$n==1)
# m1 = length(ind)
# N_k = subset$n
# N_k = sort(N_k)
# k = length(N_k)
# m1_bar <- round(runif(1,0,m1),digits=0)
#sigma <- 0.17
#psi <- log(sigma)-log(1-sigma)
#theta <- 25
#lambda <- log(theta)

library(dplyr)

load("par.RData")
n = par$n
k = par$k
N_k = par$N_k

m1 = par$m1
m1_bar = par$m1_bar

sigma = par$sigma
theta = par$theta
beta = par$beta

a <- list(theta = 2, sigma = 1, beta=1)
b <- list(theta = 0.02, sigma = 1, beta=1)

######## FUNZIONI SIGMA ########

f_psi = function(lambda,psi,m1_bar){
  S=0
  for(j in 1:(k-m1)){
    temp = lgamma(N_k[j] - (exp(psi)/(1 + exp(psi)))) - lgamma(1/(1+exp(psi)))
    S = S + temp
  }
  
  val = dbeta((exp(psi)/(1+exp(psi))),a$sigma,b$sigma, log=TRUE) +
    ((k - m1_bar)*(psi - log(1 + exp(psi)))) + S + 
    lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) -
    lgamma((exp(lambda)/exp(psi)) + exp(lambda))
  
  return(val)  
}

mh_psi = function(R,alpha,lambda,psi,m1_bar){
  psi_vec = numeric(R)
  accept=0
  
  for(i in 1:R){
    
    psi_s = rnorm(1,psi,alpha)
    
    A = min(1,exp( f_psi(lambda,psi_s,m1_bar)   
                   - f_psi(lambda,psi,m1_bar) ))
    
    if(runif(1) < A){
      psi = psi_s
      accept = accept + 1
    }
    psi_vec[i] = psi
  }
  
  return(list(a=mean(psi_vec),b=accept/R))
  
}

######## FUNZIONI THETA #########

g_lambda = function(lambda,psi,m1_bar){
  
  val = dgamma(exp(lambda),a$theta,rate=b$theta,log=TRUE) + 
    lgamma(exp(lambda)) +
    lgamma((exp(lambda)/exp(psi)) + exp(lambda) + k - m1_bar) -
    lgamma((exp(lambda)/exp(psi)) + exp(lambda)) -
    lgamma(exp(lambda) + n - m1_bar)
  
  return(val)
}

mh_lambda = function(R,alpha,lambda,psi,m1_bar){
  lambda_vec = numeric(R)
  accept = 0
  
  for(i in 1:R){
    
    lambda_s = rnorm(1,lambda,alpha)
    
    A = min(1,exp(  g_lambda(lambda_s,psi,m1_bar)  
                    - g_lambda(lambda,psi,m1_bar) ))
    
    if(runif(1) < A){
      lambda = lambda_s
      accept = accept + 1
    }
    lambda_vec[i] = lambda
  }
  
  return(list(a=mean(lambda_vec),b=accept/R))
  
}

R <- 350
it <- 500
alpha <- list(theta = 1.5, sigma = 0.2)

psi_vec_2 <- numeric(R)
sigma_accept <- numeric(R)
lambda_vec_2 <- numeric(R)
theta_accept <- numeric(R)
beta_vec <- numeric(R)
# m1_bar_vec <- numeric(R)

psi_old <- log(sigma) -  log(1-sigma)
lambda_old <- log(theta)
beta_old <- beta
m1_bar_old <- m1_bar


for (i in 1:R){
  
  psi_new = mh_psi(it,alpha$sigma,lambda_old,psi_old,m1_bar_old)$a
  sigma_accept[i] = mh_psi(it,alpha$sigma,lambda_old,psi_old,m1_bar_old)$b
  
  lambda_new = mh_lambda(it,alpha$theta,lambda_old,psi_new,m1_bar_old)$a
  theta_accept[i] = mh_lambda(it,alpha$theta,lambda_old,psi_new,m1_bar_old)$b
  
  beta_new = rbeta(1, a$beta + n - m1_bar_old, b$beta + m1_bar_old)
  
#  m1_bar_new = exp(lchoose(m1,m1_bar_old) + (n-m1_bar_old)*log(beta_new) + (m1_bar_old)*(1-beta_new)) *
#               (exp(psi_new)/(1+exp(psi_new)))^(k-m1_bar_old) * 
#               exp(lgamma((exp(lambda_new)/exp(psi_new)) + exp(lambda_new) + k - m1_bar_old) -
#               lgamma(exp(lambda_new) + n - m1_bar_old))
  
  
  psi_vec_2[i] = psi_new
  lambda_vec_2[i] = lambda_new
  beta_vec[i] = beta_new
  # m1_bar_vec[i] = m1_bar_new
  
  psi_old = psi_new
  lambda_old = lambda_new
  beta_old = beta_new
  # m1_bar_old = round(m1_bar_new,digits=0)
  
}


hist(psi_vec_2,100)
hist(lambda_vec_2,100)
hist(beta_vec,100)
# hist(m1_bar_vec,100)

plot(psi_vec_2,type="l")
plot(lambda_vec_2,type="l")
plot(beta_vec,type='l')
# plot(m1_bar_vec, type='l')

psi = mean(psi_vec_2)
lambda = mean(lambda_vec_2)

sigma = exp(psi)/(1+exp(psi))
theta = exp(lambda)
beta = mean(beta_vec)
# m1_bar = mean(m1_bar_vec)

param = list(sigma=sigma,theta=theta,beta=beta,m1_bar=m1_bar)
param
