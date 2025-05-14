#
#path = dirname(rstudioapi::getSourceEditorContext()$path)
source("NExON/functions.R")
source("NExON/NExON_func.R")

#
Tn <- 5
ts <- 1:Tn
P <- 20

# at the first time point
#
set.seed(1234)
A0 <- matrix(rbinom(P*P, 1, 0.03), nrow = P)
A0[lower.tri(A0)] = t(A0)[lower.tri(A0)]
diag(A0) <- 0
print(sum(A0[upper.tri(A0)])/sum(upper.tri(A0)))
print(A0)

#
# 3 edges become present over time; 1 becomes absent over time
#
edge_id <- which(A0 == T & upper.tri(A0), arr.ind = T)
no_edge_id <- which(A0 == F & upper.tri(A0), arr.ind = T)
npos <- 3; nneg <- 1
pos_id <- matrix(no_edge_id[sample(1:nrow(no_edge_id),npos),], ncol=2)
neg_id <- matrix(edge_id[sample(1:nrow(edge_id),nneg),], ncol=2)
pos_id_sign_change <- (ceiling(Tn/2)-1):(ceiling(Tn/2)+1)
neg_id_sign_change <- ceiling(Tn/2)
As <- replicate(Tn, A0, simplify=FALSE)

for (i in 1:npos) {
   for (t in (pos_id_sign_change[i] + 1):Tn) {
    As[[t]][pos_id[i,1],pos_id[i,2]] <-  As[[t]][pos_id[i,2],pos_id[i,1]] <- T
  }
}

for (i in 1:nneg) {
  for (t in (neg_id_sign_change[i] + 1):Tn) {
    As[[t]][neg_id[i,1],neg_id[i,2]] <- As[[t]][neg_id[i,2],neg_id[i,1]] <- F
  }
}
sparsity <- sapply(As, function(x)sum(x[upper.tri(x)])/sum(upper.tri(x)))
print(sparsity)


# simulate data
#
Ns <- sample(c(100,150,200), Tn, replace = T)
ans <- lapply(1:Tn, function(t)
  generate_network_tan(n = Ns[t],
                       p = P,
                       A = As[[t]], 
                       empirical = T))
Ys <- lapply(ans, function(x)x$Y)

# inference
#
out <- NExON(Ys,
           ts,
           debug = T,
           v0_list = rep(0.1,5))

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

