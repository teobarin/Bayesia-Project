library(dplyr)
library(readr)

# - freq.RData: file containing the frequencies of the clusters of our dataset
# - functions.RData: file containing all functions that we use to estimate parameters
# - sigma_par: list containing the parameters of the prior distribution of sigma
# - theta_par: list containing the parameters of the prior distribution of theta
# - beta_par: list containing the parameters of the prior distribution of beta
# - alpha: list containing the variance of the proposal distribution that we use
#   in the Metropolis-Hastings to update sigma and theta
# - iters: numbers of iteration of our sampling

load("C:/.../freq.RData")
load("C:/.../functions.RData")

sigma_par = list(a=1,b=1)
theta_par = list(a=2,b=0.02)
beta_par = list(a=1,b=1)
alpha = list(psi = 0.1, lambda = 0.1)
iters = 10000

d = cPY_freq(freq,sigma_par,theta_par,beta_par,alpha,iters)
d_nc = PY_freq(freq,sigma_par,theta_par,alpha,iters)
