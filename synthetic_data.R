#
#
#### DATI SINTETICI
#
#

n <- 50000
X <- numeric(n)
beta <- 0.9
sigma <- 0.2
theta <- 100
J <- numeric(n)
X[1] <- 1
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
      t <- rdunif(1, i, 1) # genera un indice tra quelli disponibili
      while (J[t]==0){
        t <- rdunif(1, i, 1)
      }
      
      # 3: non contaminante, freq 1
      if ((COUNT == 1) ){ # se quel valore l'ho già incontrato solo una volta
        X[i] <- X[t]
        m1 <- m1-1
        J[i] <- 1
      }

      # 4: non contaminante, freq > 1
      if (COUNT > 1){ # se quel valore l'ho già incontrato più volte
        X[i] <- X[t]
        J[i] <- 1
      }

    }
  }
  
  # taglio di J, ovvero
  # controlla in X tutti quelli che hanno frequenza >1 e elimina i valori di J in tutti gli
  # indici dove compare quel numero
    
}


  
  
  
  
  
  
  
  
  
  
save(X, file = "vec.RData")






























