get_annealing_ladder_ <- function(anneal, verbose) {
  # ladder set following:
  # Importance Tempering, Robert B. Gramacy & Richard J. Samworth, pp.9-10, arxiv v4
  
  k_m <- 1 / anneal[2]
  m <- anneal[3]
  
  if (anneal[1] == 1) {
    type <- "geometric"
    
    delta_k <- k_m ^ (1 / (1 - m)) - 1
    
    ladder <- (1 + delta_k) ^ (1 - m:1)
    
  } else if (anneal[1] == 2) {
    # harmonic spacing
    
    type <- "harmonic"
    
    delta_k <- (1 / k_m - 1) / (m - 1)
    
    ladder <- 1 / (1 + delta_k * (m:1 - 1))
    
  } else if (anneal[1] == 3) {
    # linear spacing
    
    type <- "linear"
    
    delta_k <- (1 - k_m) / (m - 1)
    
    ladder <- k_m + delta_k * (1:m - 1)
  } else {
    type <- "fixed"
    ladder <- k_m
  }
  
  if (verbose != 0)
    cat(paste0("** Annealing with ", type, " spacing ** \n\n"))
  
  ladder
  
}


create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
  
}

log_sum_exp <- function(x) {
  # Computes log(sum(exp(x))
  
  if (max(abs(x)) > max(x)) {
    offset <- min(x)
  } else {
    offset <- max(x)
  }
  
  log(sum(exp(x - offset))) + offset
  
}

effective_sum <- function(x) {
  # effective sum using log_sum_exp
  xp <- x[x > 0]
  xn <- x[x < 0]
  
  if (length(xp) != 0 & length(xn) != 0) {
    
    fexp(log_sum_exp(log(xp))) - exp(log_sum_exp(log(-xn)))
    
  } else if (length(xp) == 0 & length(xn) != 0) {
    
    -exp(log_sum_exp(log(-xn)))
    
  } else if (length(xp) != 0 & length(xn) == 0) {
    
    exp(log_sum_exp(log(xp)))
    
  } else if (length(xp) == 0 & length(xn) == 0) {
    
    0
    
  }
}
