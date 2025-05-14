#####################

## delta's updates ##

#####################

update_m_deltas <- function(Omegas, m1_alphas, v0s, v1, Tn, c=1) {
  
  lapply(1:Tn, function(t){
    1 / (1 + exp(
      c * log(v1 / v0s[[t]]) +
        c  * Omegas[[t]] ^ 2 * (1 / v1 ^ 2 - 1 / v0s[[t]] ^ 2) / 2 +
        pnorm(sqrt(c) * m1_alphas[[t]], log.p = T, lower.tail = F) -
        pnorm(sqrt(c) * m1_alphas[[t]], log.p = T, lower.tail = T)
    ))
  })
}

get_dstars <- function(m_deltas, v0_list, v1) {
  
  Map(function(m_deltas, v0) {
    m_deltas / v1^2 + (1 - m_deltas) / v0^2
  }, m_deltas, v0_list)
  
}


#####################

## z's updates ##

#####################

get_m_zs <- function(m_deltas, m1_alphas, Tn, c = 1) {
  
  sqrt_c <- sqrt(c)
  lapply(1:Tn, function(t){
    log_pnorm <- pnorm(sqrt_c * m1_alphas[[t]], log.p = TRUE)
    log_1_pnorm <-
      pnorm(sqrt_c * m1_alphas[[t]], log.p = TRUE, lower.tail = FALSE)
    
    imr0 <-
      inv_mills_ratio_(0, sqrt_c * m1_alphas[[t]], log_1_pnorm, log_pnorm)
    imr1 <-
      inv_mills_ratio_(1, sqrt_c * m1_alphas[[t]], log_1_pnorm, log_pnorm)
    
    m_z <- m1_alphas[[t]]  + imr0 / sqrt_c + m_deltas[[t]] * (imr1 - imr0) / sqrt_c
    
    diag(m_z) <- 0
    
    return(m_z)
  })
  
}


get_m2_zs <- function(m_zs, m1_alphas, Tn, cst = 1) {
  
  lapply(1:Tn, function(t){
    m1_alphas[[t]] * m_zs[[t]] + 1 / cst
  })
  
}

########################

## zeta ij 's updates ##

########################

update_sig2_inv_zetaij <- function(t02, P, Tn, c = 1) {
  
   matrix(c * (Tn + 1 / t02), nrow = P, ncol = P)
  
}

update_mu_zetaij <-
  function(mu_beta, sig2_inv_zeta, m_zs, ts, n0, t02, P, c = 1) {

    c * (Reduce("+", m_zs) - mu_beta * sum(ts) + n0 / t02) / sig2_inv_zeta
    # diagonal is not useful here
  }


get_m2_zeta <- function(mu_zeta, sig2_inv_zeta) {
  
  mu_zeta ^ 2 + sig2_inv_zeta ^ (-1)
  
}

#####################

## sigma's updates ##

#####################

update_alpha_sigma <- function(a_sigma, P, c = 1) {
  
  c * (P * (P-1) / 4 + a_sigma -1) + 1
  
}

update_beta_sigma <- function(m2_beta, b_sigma, c = 1) {
  
  bool_up <- upper.tri(m2_beta)
  
  c * (sum(m2_beta[bool_up]) / 2 + b_sigma)
  
}

get_m_sig2_inv <- function(alpha_sigma, beta_sigma){
  
  alpha_sigma / beta_sigma
  
}

get_m_log_sig2_inv <- function(alpha_sigma, beta_sigma) {
  
  eps <- .Machine$double.eps
  
  digamma(alpha_sigma + eps) - log(beta_sigma + eps)
  
}

####################

## beta's updates ##

####################

update_sig2_inv_beta <- function(m_sig2_inv, ts, P, c = 1) {
  
  ans <- c *  ( m_sig2_inv + sum( ts ^ 2 ) )
  return(matrix(ans, P, P))
  
}


update_mu_beta <-
  function(sig2_inv_beta,
           m_zs,
           mu_zeta,
           ts,
           Tn,
           c = 1) {
    c * (
      Reduce("+", lapply(1:Tn, function(t) ts[t] * m_zs[[t]])) -
        mu_zeta * sum(ts)
    ) / sig2_inv_beta
    
  }

get_m2_beta <- function(mu_beta, sig2_inv_beta) {
  
  (mu_beta ^ 2 + sig2_inv_beta ^ (-1)) 
  
}


####################

## omega's updates ##

####################

get_omegas <- function(dstars, Ss, Omegas, lambda, Ns, P, Tn) {
  
  lapply(1:Tn, function(t){
    for (j in 1:P) {
      IOmega_nj_nj <-
        solve(Omegas[[t]][-j, -j, drop = FALSE]) 
      # implement update based on Sigma to avoid inverting here.
      
      s_j_j <- Ss[[t]][j, j]
      
      Omegas[[t]][-j, j] <-
        Omegas[[t]][j, -j] <-
        - solve((s_j_j + lambda) * IOmega_nj_nj + 
                  diag(dstars[[t]][-j, j]), Ss[[t]][-j, j])
      Omegas[[t]][j, j] <-
        Omegas[[t]][j, -j, drop = FALSE] %*% IOmega_nj_nj %*% Omegas[[t]][-j, j] +
        Ns[t] / (lambda + s_j_j)
    }
    
    return(Omegas[[t]])
  }
  )
}

#####################

## alpha's updates ##

#####################

get_m1_alphas <- function(mu_zeta, mu_beta, ts){
  
  lapply(ts, function(t) mu_zeta + mu_beta * t)
  
}

get_var_alphas <- function(sig2_inv_zeta, sig2_inv_beta, ts){
  # minus var alphas
  lapply(ts, function(t) sig2_inv_zeta^(-1) + t^2 * sig2_inv_beta^(-1))
  
}

#####################

## E log p(y|rest) ##

#####################

e_y <- function(Omegas, Ss, Ns, Tn){
  
  sum(sapply(1:Tn,
             function(t){
               Ns[t] * determinant(Omegas[[t]], logarithm = TRUE)$modulus[1] / 2 -
                 sum(Ss[[t]] * Omegas[[t]]) / 2
             }))
  
}

############################################

## E log p(omegas | rest) - E log q(omegas) ##

############################################

e_omegas <- function(Omegas,
                     m_deltas,
                     dstars,
                     lambda,
                     v0s,
                     v1,
                     P,
                     Tn){
  
  bool_up <- upper.tri(Omegas[[1]])
  
  sum(sapply(1:Tn, function(t){
    - lambda * sum(diag(Omegas[[t]])) / 2 -
      log(v1) * sum(m_deltas[[t]][bool_up]) -
      log(v0s[[t]]) * sum(1-m_deltas[[t]][bool_up]) -
       sum(Omegas[[t]][bool_up]^2 * dstars[[t]][bool_up])/2
  }))
}

######################################################

## E log p(deltas, zs | rest) - E log q(deltas, zs) ##

######################################################
#' @importFrom stats pnorm
#' 
e_deltas_zs <- function(m_deltas,
                        m1_alphas,
                        var_alphas,
                        Tn,
                        c){
  
  eps <- .Machine$double.eps
  bool_up <- upper.tri(m_deltas[[1]])
  
  sum(sapply(1:Tn, function(t){
    var_alphas[[t]][bool_up]/2 + 
      m_deltas[[t]][bool_up] * stats::pnorm(sqrt(c) * m1_alphas[[t]][bool_up], log.p = T, lower.tail = T) + 
      (1-m_deltas[[t]][bool_up]) * stats::pnorm(sqrt(c) * m1_alphas[[t]][bool_up], log.p = T, lower.tail = F) -
      m_deltas[[t]][bool_up] * log(m_deltas[[t]][bool_up] + eps) -
      (1-m_deltas[[t]][bool_up]) * log(1-m_deltas[[t]][bool_up] + eps)
  }
  )
  )
  
}


##########################################

## E log p(zeta | rest) - E log q(zeta) ##

##########################################

e_zetaij <- function(mu_zeta,
                     sig2_inv_zeta,
                     m2_zeta,
                     n0,
                     t02){
  
  eps <- .Machine$double.eps
  bool_up <- upper.tri(mu_zeta)
  - sum(m2_zeta[bool_up])/(2 * t02) +
    n0 * sum(mu_zeta[bool_up])/t02 +
    sum((1 + log(2 * pi * sig2_inv_zeta[bool_up]^(-1) + eps))/2)
  
}

##########################################

## E log p(beta | rest) - E log q(beta) ##

##########################################

e_beta <- function(m_sig2_inv,
                   m_log_sig2_inv,
                   m2_beta,
                   sig2_inv_beta,
                   P
){
  
  eps <- .Machine$double.eps
  bool_up <- upper.tri(m2_beta)
  P * (P-1) * m_log_sig2_inv /4 -
    sum(m2_beta[bool_up])/2 * m_sig2_inv +
    sum((1 + log(2 * pi * sig2_inv_beta[bool_up]^(-1) + eps))/2 ) 
  
}


############################################

## E log p(sigma | rest) - E log q(sigma) ##

############################################

e_sigma <- function(alpha_sigma,
                    beta_sigma,
                    m_sig2_inv,
                    m_log_sig2_inv,
                    a_sigma,
                    b_sigma){
  
  eps <- .Machine$double.eps
  
  (a_sigma - alpha_sigma) * m_log_sig2_inv -
    (b_sigma - beta_sigma) * m_sig2_inv -
    alpha_sigma * log(beta_sigma + eps) + lgamma(alpha_sigma)
  
}


##########

## ELBO ##

##########

get_elbo_zetaij <- function(Omegas,
                            m_deltas,
                            mu_zeta,
                            sig2_inv_zeta,
                            m2_zeta,
                            alpha_sigma,
                            beta_sigma,
                            m_log_sig2_inv,
                            m_sig2_inv,
                            m2_beta,
                            sig2_inv_beta,
                            m1_alphas,
                            var_alphas,
                            dstars,
                            Ss,
                            lambda,
                            v0s,
                            v1,
                            n0,
                            t02,
                            a_sigma,
                            b_sigma,
                            Ns,
                            P,
                            Tn,
                            c = 1) {
  
  c * (
    e_y(Omegas, Ss, Ns, Tn) +
      e_omegas(Omegas,
               m_deltas,
               dstars,
               lambda,
               v0s,
               v1,
               P,
               Tn) +
      e_deltas_zs(m_deltas,
                  m1_alphas,
                  var_alphas,
                  Tn,
                  c) +
      e_zetaij(mu_zeta,
               sig2_inv_zeta,
               m2_zeta,
               n0,
               t02) +
      e_beta(m_sig2_inv,
             m_log_sig2_inv,
             m2_beta,
             sig2_inv_beta,
             P) +
      e_sigma(alpha_sigma,
              beta_sigma,
              m_sig2_inv,
              m_log_sig2_inv,
              a_sigma,
              b_sigma)
  )
}
# check the presence of character name in the list x; if not, set to default.
#
set_default <- function(x, name, default){
  if (!name %in% names(x) | (name %in% names(x) & is.null(x[[name]]))) {
    x[[name]] <- default
  }
  return(x)
}

inv_mills_ratio_ <- function(delta, U, log_1_pnorm_U, log_pnorm_U) {
  stopifnot(delta %in% c(0, 1))
  
  # writing explicitly the formula for pnorm(, log = TRUE) is faster...
  if (delta == 1) {
    m <- exp(-U ^ 2 / 2 - log(sqrt(2 * pi)) - log_pnorm_U)
    
    m[m < -U] <- -U[m < -U] # to do correct in other packages
    
  } else {
    m <- -exp(-U ^ 2 / 2 - log(sqrt(2 * pi)) - log_1_pnorm_U)
    
    m[m > -U] <- -U[m > -U]
    
  }
  
  m
  
}

#' @importFrom stats setNames
create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
  
}




generate_network_tan <- function(n,
                                 p,
                                 prob_edge,
                                 prob_edge_hubs,
                                 n_hubs,
                                 vec_magnitude = c(0.25, 0.75),
                                 n_connected_comp = 1,
                                 bool_scale = TRUE,
                                 empirical = F,
                                 A = NULL) {
  # same generation procedure as in Learning Graphical Models With Hubs, JMLR, 2014, p.3307
  # might be useful to compare to the implementation of the function HubNetwork from the R package hglasso (similar? same authors?)
  
  ind_hubs <- NA
  if (is.null(A)) {
    stopifnot(n_hubs <= p)
    stopifnot(n_connected_comp <= p)
    stopifnot(length(prob_edge) == 1)
    stopifnot(length(prob_edge_hubs) %in% c(1, n_hubs)) # probability of edge will be hub-specific if a vector of length n_hubs is supplied
    if (length(prob_edge_hubs) == 1) {
      prob_edge_hubs <- rep(prob_edge_hubs, n_hubs)
    }
    
    check_zero_one_(prob_edge)
    check_zero_one_(prob_edge_hubs)
    stopifnot(prob_edge < prob_edge_hubs)
    stopifnot(n_connected_comp == 1) # not implemented yet
    
    ind_hubs <- sort(sample(1:p, n_hubs))
    
    # generate adjacency matrix
    #
    A <- matrix(sample(
      c(0, 1),
      size = p ^ 2,
      prob = c(1 - prob_edge, prob_edge),
      replace = TRUE
    ),
    nrow = p,
    ncol = p)
    
    for (h in 1:n_hubs) {
      ind_h <- ind_hubs[h]
      prob_edge_h <- prob_edge_hubs[h]
      A[ind_h, ] <- sample(
        c(0, 1),
        size = p,
        prob = c(1 - prob_edge_h, prob_edge_h),
        replace = TRUE
      )
      
      A[, ind_h] <- sample(
        c(0, 1),
        size = p,
        prob = c(1 - prob_edge_h, prob_edge_h),
        replace = TRUE
      )
      
    }
    
    A[lower.tri(A)] <-
      t(A)[lower.tri(A)]  # we actually only use the upper part of A, as the adjacency matrix is symmetric
    diag(A) <- 0
  }
  
  nb_edges <- sum(A == 1)
  
  # matrix E
  E <- A
  
  E[A == 1] <- runif(nb_edges, min = vec_magnitude[1], max = vec_magnitude[2])
  
  E_bar <- (E + t(E)) / 2
  
  msign <- matrix(1, nrow = nrow(E), ncol = ncol(E))
  msign[upper.tri(msign)] <- sample(c(-1,1), size = sum(upper.tri(msign)),  prob = c(0.5, 0.5), replace = TRUE)
  msign[lower.tri(msign)] <- t(msign)[lower.tri(msign)]
  E_bar <- E_bar * msign
  
  # minimum eigenvalue
  min_eigen <- min(eigen(E_bar, only.values = TRUE)$values)
  if (min_eigen < 0) {
    Omega <- E_bar + (0.1 - min_eigen) * diag(p)
  } else{
    Omega <- E_bar + 0.1 * diag(p)
  }
  
  if(empirical){
    Y <- MASS::mvrnorm(n, rep(0, p), solve(Omega), empirical = T)
  }else{
    Y <- mvtnorm::rmvnorm(n, rep(0, p), solve(Omega))
  }
  
  
  if (bool_scale) {
    Y <- scale(Y)
  }
  
  create_named_list_(A, Omega, Y, ind_hubs, nb_edges)
}


#' @importFrom stats uniroot
get_n0_t02 <- function(p, p_star) {
  
  E_p_t <- p_star[1]
  # V_p_t <- min(p_star[2], floor(2 * p / 3))
  V_p_t <- p_star[2]
  
  dn <- 1e-6
  up <- 1e5
  
  # Get n0 and t02
  #
  tryCatch(t02 <- stats::uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of edges when no hubs appear. ",
                  "Please change their values."))
    })
  
  #
  n0 <- get_mu(E_p_t, t02, p)
  
  create_named_list_(n0, t02)
}