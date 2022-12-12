#
#
#### DATI SINTETICI
#
#
library(dplyr)

n <- 100000
X <- numeric(n)
J <- numeric(n)

beta <- 0.95
sigma <- 0.2
theta <- 100

clust <- 0
m1 <- 0
m1_bar <- 0
m1_tilde <- 0

a <- numeric(n)
b <- numeric(n)
c <- numeric(n)


for (i in 1:n){
  
  m1 <- m1_bar + m1_tilde
  
  a[i] <- (theta + (clust - m1_bar)*sigma)/(theta + i - m1_bar)
  b[i] <- (m1_tilde*(1 - sigma))/(theta + i - m1_bar)
  c[i] <- ((clust - m1)*sigma + n - m1_bar)/(theta + i - m1_bar)
  
  u <- runif(1,0,1)
  
  if (u > beta){
    # 1: specie contaminante nuova
    X[i] <- clust + 1
    J[i] <- 0
    clust <- clust + 1
    m1_bar <- m1_bar + 1
  }
  else{
    v <- runif(1,0,1)
    
    if (v < a[i]){
      # 2: specie non contaminante nuova
      X[i] <- clust + 1
      J[i] <- 1
     clust <- clust + 1
      m1_tilde <- m1_tilde + 1
    }
    else{
     
      if (v < (a[i] + b[i])){
        # 3: non contaminante, freq 1
        t <- round(runif(1, 1, clust),digits=0) # genera un indice tra quelli disponibili
        while (J[t] == 0 && length(X[t])>1){
          t <- round(runif(1, 1, clust),digits=0)
        }
        
        X[i] <- X[t]
        m1_tilde <- m1_tilde - 1
        J[i] <- 1
      }
      else{
        # 4: non contaminante, freq > 1
        t <- round(runif(1, 1, clust),digits=0) # genera un indice tra quelli disponibili
        while (J[t] == 0 && length(X[t])==1){
          t <- round(runif(1, 1, clust),digits=0)
        }
          X[i] <- X[t]
          J[i] <- 1
        }
      }
      
    }
  }
  

# lego i gruppi a J
data = as.data.frame(cbind(X,J))

# rimuovo dal dataset i singleton contaminanti
ind = which(data$J==0)
data_no_cont = data[-ind,]

# trovo la numerosità di ogni gruppo
subset = count(data_no_cont,X)

# Rimuovo i singleton non cantaminanti dal dataset e genero il vettore con le numerosità dei gruppi
ind_2 = which(subset$n==1)
subset_2 = subset[-ind_2,]
N_k = subset_2$n
N_k = sort(N_k)

# numero di gruppi
k = m1 + length(N_k)

# salvo tutti i parametri che mi servono per controllare se il modello funziona
par = list(n=n, m1=m1, m1_bar=m1_bar, N_k=N_k, k=k, beta=beta, sigma=sigma, theta=theta)
par

save(par, file = "par.RData")
