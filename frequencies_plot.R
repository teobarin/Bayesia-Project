library(dplyr)
library(reader)

data <- read.csv("india_small.csv", sep ="\t", header=T)
freq <- as.vector(table(data$verbatimScientificName))

N_k = freq 
N_k = sort(N_k)
k = length(N_k)

fr = as.data.frame(table(as.data.frame(N_k)))
plot(log(as.numeric(fr$Var1)),log(fr$Freq),xlab="n_j(log-scale)",ylab="m_j(log-scale)")
