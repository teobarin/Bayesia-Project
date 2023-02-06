library(dplyr)
library(readr)

load("C:/Users/marco/Desktop/Univerisità/Bayesian Statistics/Progetto/freq.RData")
load("C:/Users/marco/Desktop/Univerisità/Bayesian Statistics/Progetto/functions.RData")

sigma_par = list(a=1,b=1)
theta_par = list(a=2,b=0.02)
beta_par = list(a=1,b=1)
alpha = list(psi = 0.1, lambda = 0.1)
iters = 10000

d = cPY_freq(freq,sigma_par,theta_par,beta_par,alpha,iters)
d_nc = PY_freq(freq,sigma_par,theta_par,alpha,iters)
