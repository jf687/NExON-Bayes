update_m_delta <- function(Omega, m_tau, m_log_rho, m_log_one_minus_rho, v0, v1, c=1) {
  1 / (1 + exp(
    c * log(v1 / v0) +
      c * m_tau * Omega ^ 2 * (1 / v1 ^ 2 - 1 / v0 ^ 2) / 2 +
      c * m_log_one_minus_rho -
      c * m_log_rho
  ))
}

# v0
# ===
update_a_v0 <- function(a, m_delta){
  bool_up <- upper.tri(m_delta)
  0.5*effective_sum(1 - m_delta[bool_up]) + a
}

update_b_v0 <- function(b, m_delta, Omega){
  bool_up <- upper.tri(m_delta)
  0.5*sum((1 - m_delta[bool_up])*(Omega[bool_up]^2)) + b
}

get_m_v0 <- function(a_v0, b_v0) {
  bool_log_sum_exp <- F
  if (!bool_log_sum_exp) {
    b_v0 / a_v0
  } else{
    eps <- .Machine$double.eps
    exp(log(a_v0 + eps) - log(b_v0 + eps))
  }
}



# tau
# === 
update_a_tau <- function(a, p, c = 1) {
  c * (p * (p - 1) / 4 + a  - 1) + 1
}

update_b_tau <- function(Omega, E1, b, c = 1) {
  bool_up <- upper.tri(Omega)
  c * (b + sum(Omega[bool_up] ^ 2 * E1[bool_up]) / 2)
}

get_m_tau <- function(a_tau, b_tau) {
  bool_log_sum_exp <- F
  if (!bool_log_sum_exp) {
    1 #a_tau / b_tau
  } else{
    eps <- .Machine$double.eps
    exp(log(a_tau + eps) - log(b_tau + eps))
  }
}

get_m_log_tau <- function(a_tau, b_tau) {
  #eps <- .Machine$double.eps
  0
  #digamma(a_tau + eps) - log(b_tau + eps)
}


get_E1 <- function(P, v0, v1) {
  bool_log_sum_exp <- F
  # logsumexp not suitable for extremes
  if (!bool_log_sum_exp) {
    ans <- (1 - P) / v0 ^ 2 + P / v1 ^ 2
    return(ans)
  } else{
    tmp <-
      matrix(mapply(function(x, y)
        effective_sum(c(x, y)), (1 - P) / v0 ^ 2, P / v1 ^ 2), ncol = ncol(P))
    return(tmp)
  }
  # hist(ans - tmp)
}


# rho
# === 
update_a_rho <- function(m_delta, ar, c = 1) {
  # c * (m_delta + ar - 1) + 1
  bool_up <- upper.tri(m_delta)
  c * (effective_sum(m_delta[bool_up]) + ar - 1) + 1
}

update_b_rho <- function(m_delta, br, c = 1) {
  bool_up <- upper.tri(m_delta)
  c * (effective_sum(1 - m_delta[bool_up]) + br - 1) + 1
}

get_m_log_rho <- function(a_rho, b_rho) {
  eps <- .Machine$double.eps
  digamma(a_rho + eps) - digamma(a_rho + b_rho + eps)
}

get_m_log_one_minus_rho <- function(a_rho, b_rho) {
  eps <- .Machine$double.eps
  digamma(b_rho + eps) - digamma(a_rho + b_rho + eps)
}


# M-steps:
# ====== #


get_omega <- function(E1, S, Omega, lambda, n, p) {
  for (j in 1:p) {
    IOmega_nj_nj <-
      solve(Omega[-j, -j, drop = FALSE]) # implement update based on Sigma to avoid inverting here.
    
    s_j_j <- S[j, j]
    
    Omega[-j, j] <-
      Omega[j, -j] <-
      -solve((s_j_j + lambda) * IOmega_nj_nj + diag(E1[-j, j]), S[-j, j])
    Omega[j, j] <-
      Omega[j, -j, drop = FALSE] %*% IOmega_nj_nj %*% Omega[-j, j] + n / (lambda + s_j_j)
    
  }
  
  Omega
  
}

get_elbo_GM <- function(Omega,
                                  m_delta,
                                  a_tau,
                                  b_tau,
                                  m_tau,
                                  m_log_tau,
                                  a_rho,
                                  b_rho,
                                  m_log_rho,
                                  m_log_one_minus_rho,
                                  S,
                                  E1,
                                  lambda,
                                  v0,
                                  v1,
                                  a,
                                  b,
                                  ar,
                                  br,
                                  n,
                                  p) {
  eps <- .Machine$double.eps
  bool_up <- upper.tri(Omega)
  #
  n * determinant(Omega, logarithm = TRUE)$modulus[1] / 2 -
    sum(S * Omega) / 2 -
    lambda * sum(diag(Omega)) / 2 -
    log(v1) * sum(m_delta[bool_up]) -
    log(v0) * sum(1-m_delta[bool_up]) -
    m_tau * (sum(Omega[bool_up]^2 * E1[bool_up])/2 + b - b_tau) +
    m_log_tau * (p * (p-1)/4 + a - a_tau) +
    (sum(m_delta[bool_up]) + ar - a_rho) * m_log_rho +
    (sum(1 - m_delta[bool_up]) + br - b_rho) * m_log_one_minus_rho -
    sum(m_delta[bool_up] * log(m_delta[bool_up] + eps)) -
    sum((1-m_delta[bool_up]) * log(1-m_delta[bool_up] + eps)) - 
    a_tau * log(b_tau + eps) + lgamma(a_tau) +
    lgamma(a_rho) + lgamma(b_rho) - lgamma(a_rho + b_rho)
}


