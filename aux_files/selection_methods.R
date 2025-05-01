

network_selection_val <- function(t, Omega, m_delta, Y, gamma = 0.15, method = "eBIC") {
  
  N <- nrow(Y)
  
  Omega[m_delta <= t] <- 0
  
  d <- sum(Omega[lower.tri(Omega)] != 0)
  
  det.Omega <- determinant(Omega, logarithm = TRUE)$modulus[1]
  if(method == "AIC"){
    return( N*(sum(diag(cov(Y) %*% Omega)) -  det.Omega) + 2*d)
  }
  
  if(method == "BIC"){
    return(N*(sum(diag(cov(Y) %*% Omega)) -  det.Omega) + log(N) * d)
  }
  
  if(method == "eBIC"){
    return(N*(sum(diag(cov(Y) %*% Omega)) -  det.Omega) + log(N) * d + 4*gamma*d*log(ncol(Y)))
  }
}



par.func <- function(i, v0_range, Y ,gamma){
  load_list_hyper(ncol(Y))
  v_0_ <- v0_range[i]
  out <- GM(Y, list_hyper, set_v0 = v_0_)
  return(network_selection_val(0.5, out$Omega, out$m_delta, Y, method = "eBIC", gamma = gamma))
}







select_optimal_v0_t <- function(Y, v0_range = 10^seq(log10(0.02), log10(0.2), length.out = 30),
                                     plot_ = TRUE, gamma = 0.15, save_outs = F, t_range = seq(0,1, length.out = 1000)) {
  
  list_hyper <- list(lambda = 2,
                     v0 = 0.5,
                     v1 = 100,
                     a = 2,
                     b = 2, 
                     ar = 1,
                     br = ncol(Y))
  
  outs <- list()
  IC_values <- list()
  for (i in seq_along(v0_range)){
    v_0_ <- v0_range[i]
    print(i)
    
    
    out <- GM(Y, list_hyper, set_v0 = v_0_)
    IC_values[[i]] <- list()
    for(t in seq_along(t_range)){
      IC_values[[i]][[t]] <- network_selection_val(t_range[[t]], out$Omega, out$m_delta, Y, gamma = gamma, method = "eBIC")
    }
    if(save_outs){
      outs[[i]] <- out
    }
    
  }
  IC_values <- matrix(unlist(IC_values), nrow = length(t_range))
  opt_v0_ind <- which(IC_values == min(IC_values), arr.ind = T)[[1,2]]
  opt_t_ind <- which(IC_values == min(IC_values), arr.ind = T)[[1,1]]
  
  cat("minimum BIC: ", IC_values[opt_t_ind,opt_v0_ind],"\n")
  optimal_v0 <- v0_range[opt_v0_ind]
  optimal_t <- v0_range[opt_t_ind]
  cat("optimal v0 :", optimal_v0,"\n")
  cat("optimal t :", optimal_t,"\n")
  
  
  
  
  #v_0 against eBIC values
  if(plot_){
  plot(v0_range, IC_values[50,], type = "b", log = "x", 
         xlab = expression(v[0]), ylab = "eBIC", 
         main = paste("eBIC vs",expression(v[0])),
         col = "blue", pch = 19)
    abline(v = optimal_v0, col = "red", lty = 2)
    legend("topright", legend = paste("Optimal v_0 =", round(optimal_v0, 4)),
           col = "red", lty = 2)
  }
  
  
  return(list(optimal_v0=optimal_v0, optimal_t = optimal_t, IC_values = IC_values, outs = outs))
}

