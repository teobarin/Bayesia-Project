library(dplyr)
library(ggplot2)
library(readr)

load("C:/Users/marco/Desktop/Univerisità/Bayesian Statistics/Progetto/pred_100.RData")

col <- c("dodgerblue", "darkorange2","gray48")

########### log-log ##############
n = sum(freq)
M_l <- as.numeric(table(factor(freq, levels = 1:n))) # numero di cluster con frequenza l=1,...,n
idx <- 1:(which.min(M_l)+30) # which(M_l > 0)

expected_m_py <- Vectorize(function(m, n, sigma, alpha) {
  out <- log(alpha) + lchoose(n, m) + lgamma(m - sigma) - lgamma(1 - sigma) - lgamma(alpha + n) + lgamma(alpha) + lgamma(alpha + sigma + n - m) - lgamma(alpha + sigma)
  exp(out)
}, vectorize.args = "m")

est_df_temp <- data.frame(
  x = idx, 
  y = expected_m_py(idx, n = sum(freq), sigma = quantile(d$sigma, 0.5), alpha = quantile(d$theta, 0.5)), 
  y_low = expected_m_py(idx, n = sum(freq), sigma = quantile(d$sigma, 0.05), alpha = quantile(d$theta, 0.05)), 
  y_up = expected_m_py(idx, n = sum(freq), sigma = quantile(d$sigma, 0.95), alpha = quantile(d$theta, 0.95)),
  yPY = expected_m_py(idx, n = sum(freq), sigma = quantile(d_nc$sigma_nc, 0.5), alpha = quantile(d_nc$theta_nc, 0.5)), 
  y_lowPY = expected_m_py(idx, n = sum(freq), sigma = quantile(d_nc$sigma_nc, 0.05), alpha = quantile(d_nc$theta_nc, 0.05)), 
  y_upPY = expected_m_py(idx, n = sum(freq), sigma = quantile(d_nc$sigma_nc, 0.95), alpha = quantile(d_nc$theta_nc, 0.95))
)

p <- ggplot(data.frame(x = idx[M_l != 0], y = M_l[idx[M_l != 0]]), aes(x = x, y = y)) +  
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_line(data = est_df_temp,
            mapping = aes(x = x, y = yPY), col = col[2]) +
  geom_ribbon(data = est_df_temp,
              mapping = aes(x = x, ymin = y_lowPY, ymax = y_upPY), fill = col[2], alpha = 0.2) +
  geom_line(data = est_df_temp,
            mapping = aes(x = x, y = y), col = col[1]) +
  geom_ribbon(data = est_df_temp,
              mapping = aes(x = x, ymin = y_low, ymax = y_up), fill = col[1], alpha = 0.2) +
  geom_segment(data = est_df_temp, mapping = aes(x = 1, xend = 2, y = y[1] + colMeans(cbind(d$sigma,d$theta,d$m1_bar))[3] , 
                                                 yend = y[2]), col = col[3], lty = 2) +
  geom_point(shape = 1, alpha = 0.9, col = "black") + 
  geom_point(shape = 16, alpha = 0.3, col = "black") +
  xlab(expression(n[j] * ' (log-scale)')) + 
  ylab(expression(m[j] * ' (log-scale)'))

p
################################################################################
d_sigma = data.frame(x = d$sigma, y = d_nc$sigma_nc)

p_sigma <- ggplot(d_sigma, aes(x = x)) +
  geom_density(aes(fill = "cPY"), col = col[2], alpha = 0.6,bw = 0.025) +
  geom_density(aes(x = y, fill = "PY"),col=col[1], alpha = 0.6, bw = 0.025) +
  scale_fill_manual(values = c("cPY" = col[2], "PY" = col[1]), name = "Model") +
  geom_vline(xintercept = mean(d_sigma$x), linetype="dashed", 
             color = col[3], size=1) +
  geom_vline(xintercept = mean(d_sigma$y), linetype="dashed", 
             color = col[3], size=1) +
  ggtitle(expression(sigma)) +
  xlab("Value") +
  ylab("Density") +
  theme_bw()+
  xlim(min(colMeans(d_sigma))-0.2,max(colMeans(d_sigma))+0.2)+
  theme(plot.title = element_text(size = 20,hjust=0.5))

p_sigma
################################################################################
d_theta = data.frame(x = d$theta, y = d_nc$theta_nc)

p_theta <- ggplot(d_theta, aes(x = x)) +
  geom_density(aes(fill = "cPY"), col = col[2], alpha = 0.6,bw=2.5) +
  geom_density(aes(x = y, fill = "PY"),col=col[1], alpha = 0.6,bw=2.5) +
  scale_fill_manual(values = c("cPY" = col[2], "PY" = col[1]), name = "Model") +
  geom_vline(xintercept = mean(d_theta$x), linetype="dashed", 
             color = col[3], size=1) +
  geom_vline(xintercept = mean(d_theta$y), linetype="dashed", 
             color = col[3], size=1) +
  ggtitle(expression(vartheta)) +
  xlab("Value") +
  ylab("Density") +
  theme_bw()+
  xlim(min(colMeans(d_theta))-20,max(colMeans(d_theta))+20)+
  theme(plot.title = element_text(size = 20,hjust=0.5))

p_theta
################################################################################
d_m1 = data.frame(x = m1_pred_vec, y = m1_pred_nc_vec)

p_m1 <- ggplot(d_m1, aes(x = x)) +
  geom_density(aes(fill = "cPY"), col=col[2],alpha = 0.6,bw=2.5) +
  geom_density(aes(x = y, fill = "PY"), col=col[1],alpha = 0.6,bw=2.5) +
  scale_fill_manual(values = c("cPY" = col[2], "PY" = col[1]), name = "Model") +
  geom_vline(xintercept = mean(d_m1$x), linetype="dashed", 
             color = col[3], size=1) +
  geom_vline(xintercept = mean(d_m1$y), linetype="dashed", 
             color = col[3], size=1) +
  ggtitle(expression(E*"["*N[m][","][1]^(n)*"|"*X[1]*",...,"*X[n]*","*M[n]*"]"))+
  xlab("Value") +
  ylab("Density") +
  theme_bw()+
  xlim(min(colMeans(d_m1))-10,max(colMeans(d_m1))+10)+
  theme(plot.title = element_text(size = 15,hjust=0.5))

p_m1
################################################################################
d_k = data.frame(x = k_pred_vec, y = k_pred_nc_vec)

p_k <- ggplot(d_k, aes(x = x)) +
  geom_density(aes(fill = "cPY"), col=col[2],alpha = 0.6,bw=2.5) +
  geom_density(aes(x = y, fill = "PY"), col=col[1],alpha = 0.6,bw=2.5) +
  scale_fill_manual(values = c("cPY" = col[2], "PY" = col[1]), name = "Model") +
  geom_vline(xintercept = mean(d_k$x), linetype="dashed", 
             color = col[3], size=1) +
  geom_vline(xintercept = mean(d_k$y), linetype="dashed", 
             color = col[3], size=1) +
  ggtitle(expression(E*"["*K[m]^(n)*"|"*X[1]*",...,"*X[n]*","*M[n]*"]"))+
  xlab("Value") +
  ylab("Density") +
  theme_bw()+
  xlim(min(colMeans(d_k))-10,max(colMeans(d_k))+10)+
  theme(plot.title = element_text(size = 15,hjust=0.5))

p_k

