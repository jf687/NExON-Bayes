install.packages("renv")
library(renv)

renv::restore()

source("aux_files/net_gen.R")
source("aux_files/performance_metrics.R")
source("aux_files/selection_methods.R")
source("aux_files/make_plots.R")
source("BayCON/BayCON_func.R")
source("BayCON/functions.R")
source("SSL/fun_GM.R")
source("SSL/fun_utils.R")
source("SSL/GM.R")

# Generate data from known truths
ts <- 1:4

nets <- create_Cn_networks(P = 50, Cn = 4, seed_1 = 252, seed_2 = 123, 
                              frac_change = 0.4, Ns_sample = c(150,150,150))
Ys <- nets$Ys

v0_list <- list()

# Select optimal v0s
for(c in 1:length(ts)){
  v0_list[[c]] <- select_optimal_v0_t(Ys[[c]], gamma = 0.35)$optimal_v0
}

# Run Inference Algorithm
out <- BayCON(Ys, ts, v0_list = v0_list, debug = T)
out_nets<- lapply(out$estimates$m_deltas, function(x) abs(x)> 0.5)


eval_nets <- evaluate_network_list(nets$As, out_nets,0.5)
conf.mat <- eval_nets$conf.mat
cat("Precision achieved:\n",precision(conf.mat),"\n")
cat("Recall achieved:\n",recall(conf.mat),"\n")

# Plot the estimated networks

layout <- create_layout(nets$As[[1]])
create_network_plots(out_nets, layout_coords = layout)

# Plot the true networks


create_network_plots(nets$As, layout_coords = layout,
                     color = "blue", edge_color = "green")

par(mfrow = c(1,2))
plot(out$debugs$vec_ELBO_CM, xlab = "iterations", ylab = "ELBOs at M-step")
plot(unlist(out$debugs$list_ELBO), xlab = "iterations", ylab = "ELBOs")

data <- out$estimates$mu_beta[upper.tri(out$estimates$mu_beta)]
data_frame <- data.frame(data)

ggplot(data = data_frame, aes(x = data)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "skyblue", alpha = 0.7) +
  stat_bin(binwidth = 0.05, geom = "text", aes(label = ..count..), vjust = -0.5) +
  labs(x = "(b) beta", y = "Frequency", title = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = 'gray'),
    panel.grid.minor = element_line(size = 0.25, linetype = 'dashed', color = 'gray')
  )

data <- out$estimates$mu_zeta[upper.tri(out$estimates$mu_zeta)]
data_frame <- data.frame(data)

ggplot(data = data_frame, aes(x = data)) +
  geom_histogram(binwidth = 0.5, color = "black", fill = "navyblue", alpha = 0.7) +
  stat_bin(binwidth = 0.5, geom = "text", aes(label = ..count..), vjust = -0.25) +
  labs(x = "(a) zeta", y = "Frequency", title = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = 'gray'),
    panel.grid.minor = element_line(size = 0.25, linetype = 'dashed', color = 'gray')
  )
out$estimates$mu_beta[pos_id]
out$estimates$mu_beta[neg_id]

sapply(ans, function(x)x$A[pos_id])
sapply(out$estimates$m_deltas, function(x)x[pos_id])

sapply(ans, function(x)x$A[neg_id])
sapply(out$estimates$m_deltas, function(x)x[neg_id])


cat("Precision achieved:\n",precision(conf.mat),"\n")
cat("Recall achieved:\n",recall(conf.mat),"\n")


