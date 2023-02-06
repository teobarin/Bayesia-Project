library(dplyr)
library(readr)

# - freq.RData is the file containing the frequencies of the clusters of our dataset
# - functions.RData is the file containing all functions that we use to estimate parameters
# - sigma_par is a list containing the parameters of the prior distribution of sigma
# - theta_par is a list containing the parameters of the prior distribution of theta
# - beta_par is a list containing the parameters of the prior distribution of beta
# - alpha is a list containing the variance of the proposal distribution that we use
#   in the Metropolis-Hastings to update sigma and theta
# - iters are the number of iteration of our sampling
# - cPY_freq is the function that returned the vectors containing the values 
#   of the parameters sampled with the contaminated Pitma-Yor updating functions
# - PY_freq is the function that returned the vectors containing the values 
#   of the parameters sampled with the Pitma-Yor updating functions   

load("C:/.../freq.RData")
load("C:/.../functions.RData")

sigma_par = list(a=1,b=1)
theta_par = list(a=2,b=0.02)
beta_par = list(a=1,b=1)
alpha = list(psi = 0.1, lambda = 0.1)
iters = 10000

d = cPY_freq(freq,sigma_par,theta_par,beta_par,alpha,iters)
d_nc = PY_freq(freq,sigma_par,theta_par,alpha,iters)