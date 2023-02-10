library(dplyr)
#j = 2

data <- read.csv("india_small.csv", sep ="\t", header=T)
freq <- as.vector(table(data$verbatimScientificName))

#data = read.csv(paste0("synthetic_data_", j,".csv"), header = TRUE)
#n = dim(data)[1]

#subset = count(data, data[2])[2]
N_k = freq #frequenze gruppi
N_k = sort(N_k)
k = length(N_k)

fr = as.data.frame(table(as.data.frame(N_k)))
plot(log(as.numeric(fr$Var1)),log(fr$Freq),xlab="n_j(log-scale)",ylab="m_j(log-scale)")
