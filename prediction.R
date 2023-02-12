library(dplyr)
library(readr)

load("C:/Users/marco/Desktop/Univerisità/Bayesian Statistics/Progetto/freq.RData")
load("C:/Users/marco/Desktop/Univerisità/Bayesian Statistics/Progetto/functions.RData")

sigma_par = list(a=1,b=1)
theta_par = list(a=2,b=0.02)
beta_par = list(a=1,b=1)
alpha = list(psi = 0.1, lambda = 0.1)
iters = 10000

it = 250
p = 0.2
m = round(p*sum(freq),digits=0)
data = extract_data(freq)

l = create_freq(it,data,m)

fr = l$fr
m1_true_vec = l$m1
k_true_vec = l$k

m1_pred_vec = numeric(it)
k_pred_vec = numeric(it)

m1_pred_nc_vec = numeric(it)
k_pred_nc_vec = numeric(it)

err_m1_vec = numeric(it)
err_k_vec = numeric(it)

err_m1_nc_vec = numeric(it)
err_k_nc_vec = numeric(it)

m1_true = mean(m1_true_vec)
k_true = mean(k_true_vec)


for(i in 1:it){
  
  f_tr = fr[[i]]
  
  d_tr = cPY_freq(f_tr,sigma_par,theta_par,beta_par,alpha,iters)
  d_tr_nc = PY_freq(f_tr,sigma_par,theta_par,alpha,iters)
  
  m1_pred = pred_singl(d_tr$n,m,mean(d_tr$sigma),mean(d_tr$theta),mean(d_tr$beta),d_tr$k,mean(d_tr$m1_bar))
  k_pred = pred_clust(d_tr$n,m,mean(d_tr$sigma),mean(d_tr$theta),mean(d_tr$beta),d_tr$k,mean(d_tr$m1_bar))
  m1_pred_vec[i] = m1_pred
  k_pred_vec[i] = k_pred
  
  m1_pred_nc = pred_singl_nc(d_tr_nc$n,m,mean(d_tr_nc$sigma_nc),mean(d_tr_nc$theta_nc),d_tr$k)
  k_pred_nc = pred_clust_nc(d_tr_nc$n,m,mean(d_tr_nc$sigma_nc),mean(d_tr_nc$theta_nc),d_tr$k)
  m1_pred_nc_vec[i] = m1_pred_nc
  k_pred_nc_vec[i] = k_pred_nc
  
  err_m1_vec[i] = abs(m1_true_vec[i]- m1_pred)
  err_k_vec[i] = abs(k_true_vec[i] - k_pred)
  
  err_m1_nc_vec[i] = abs(m1_true_vec[i]- m1_pred_nc)
  err_k_nc_vec[i] = abs(k_true_vec[i] - k_pred_nc)
  
}

e_m1 = mean(err_m1_vec)
e_k = mean(err_k_vec)
e_m1_nc = mean(err_m1_nc_vec)
e_k_nc = mean(err_k_nc_vec)

