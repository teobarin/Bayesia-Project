#
#
#### DATI SINTETICI
#
#
library(dplyr)

n <- 50000
X <- numeric(n)
beta <- 0.9
sigma <- 0.2
theta <- 100
J <- numeric(n)
X[1] <- 1
J[1] <- 1 #inizializzo J[1] con 1 perchè con J[1]=0 c'erano problemi 
          # nel ciclo for perchè rimaneva bloccato nel while
dummy <- 2
m1 <- 0


for (i in 2:n){
  m1_mean <- m1*(1-beta) # media di una binomiale
  u <- runif(1,0,1)
  v <- runif(1,0,1)
  
  if (u>beta){
    # 1: specie contaminante nuova
    X[i] <- dummy
    J[i] <- 0
    dummy <- dummy+1
    m1 <- m1+1
  }
  if (u <= beta){
    if (v < ((theta + (dummy-m1_mean)*sigma) / (theta+n-m1_mean))){
      # 2: specie non contaminante nuova
      X[i] <- dummy
      J[i] <- 1
      dummy <- dummy+1
      m1 <- m1+1
    }
    else{
      t <- round(runif(1, 1, i),digits=0) # genera un indice tra quelli disponibili
      while (J[t]==0){
        t <- round(runif(1, 1, i),digits=0)
        # rifalla finché non ho un valore non contaminante
      }
      
      # 3: non contaminante, freq 1
      if ((length(which(X == X[t])))==1){ # se quel valore l'ho già incontrato solo una volta
        X[i] <- X[t]
        m1 <- m1-1
        J[i] <- 1
      }
      else{
      # 4: non contaminante, freq > 1
      if (((length(which(X == X[t])))>1)){ # se quel valore l'ho già incontrato più volte
        X[i] <- X[t]
        J[i] <- 1
      }
      }
      
    }
  }
  
  # taglio di J, ovvero
  # controlla in X tutti quelli che hanno frequenza >1 e elimina i valori di J in tutti gli
  # indici dove compare quel numero
  
}

# lego i gruppi a J
data = as.data.frame(cbind(X,J))

# Trovo quanti sono i dati contaminati che rappresentano m1_bar
ind = which(data$J==0)
m1_cont = length(ind)

# rimuovo dal dataset i singleton contaminanti
data_no_cont = data[-ind,]

# trovo la numerosità di ogni gruppo
subset = count(data_no_cont,X)

# Trovo i singleton non contaminanti
ind_2 = which(subset$n==1)
m1_no_cont = length(ind_2)

# Rimuovo i singleton non cantaminanti dal dataset e genero il vettore con le numerosità dei gruppi
subset_2 = subset[-ind_2,]
N_k = subset_2$n
N_k = sort(N_k)

# numero di singleton
m1 = m1_cont + m1_no_cont

# numero di gruppi
k = m1 + length(N_k)


par = list(n=n,m1=m1,m1_bar=m1_cont,N_k=N_k,k=k)
par

save(par, file = "par.RData")
