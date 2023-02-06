library(dplyr)
library(readr)

# - it: numbers of time that we computer m1_pred and k_pred
# - p : percentage of data that we remove from dataset and we use for validation

load("C:/Users/marco/Desktop/Univerisità/Bayesian Statistics/Progetto/freq.RData")
load("C:/Users/marco/Desktop/Univerisità/Bayesian Statistics/Progetto/functions.RData")

sigma_par = list(a=1,b=1)
theta_par = list(a=2,b=0.02)
beta_par = list(a=1,b=1)
alpha = list(psi = 0.1, lambda = 0.1)
iters = 10000

it=100
p = 0.2
m = round(p*sum(freq),digits=0)
data = extract_data(freq)

l = create_freq(it,data,m)

m1_pred_vec = numeric(it)
k_pred_vec = numeric(it)

m1_pred_nc_vec = numeric(it)
k_pred_nc_vec = numeric(it)

for(i in 1:it){
  
  f_tr = l[[i]]
  
  d_tr = cPY_freq(f_tr,sigma_par,theta_par,beta_par,alpha,iters)
  d_tr_nc = PY_freq(f_tr,sigma_par,theta_par,alpha,iters)
  
  m1_pred = pred_singl(d_tr$n,m,mean(d_tr$sigma),mean(d_tr$theta),mean(d_tr$beta),d_tr$k,mean(d_tr$m1_bar))
  k_pred = pred_clust(d_tr$n,m,mean(d_tr$sigma),mean(d_tr$theta),mean(d_tr$beta),d_tr$k,mean(d_tr$m1_bar))
  m1_pred_vec[i] = m1_pred
  k_pred_vec[i] = k_pred
  
  m1_pred_nc = pred_singl(d_tr_nc$n,m,mean(d_tr_nc$sigma_nc),mean(d_tr_nc$theta_nc),mean(d_tr$beta),d_tr$k,mean(d_tr$m1_bar))
  k_pred_nc = pred_clust(d_tr_nc$n,m,mean(d_tr_nc$sigma_nc),mean(d_tr_nc$theta_nc),mean(d_tr$beta),d_tr$k,mean(d_tr$m1_bar))
  m1_pred_nc_vec[i] = m1_pred_nc
  k_pred_nc_vec[i] = k_pred_nc
}


