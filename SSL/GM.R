# Model: Graphical model with a continuous spike-and-slab prior on edge inclusion. 
# -------------------------------------------------
#' GM function.
#'
#' @examples
#' seed <- 123; set.seed(seed)
#' n <- 100
#' p <- 5
#' Y <- matrix(rnorm(n*p), nrow = n, ncol = p)
#'
#' list_hyper <- list(lambda = 2,
#'                    v0 = 0.5,
#'                    v1 = 100,
#'                    a = 2,
#'                    b = 2, 
#'                    ar = 1,
#'                    br = p)
#'                    
#' list_init <- list(a_rho = 1,
#'                   b_rho = 1,
#'                   a_tau = 1,
#'                   b_tau = 1)
#'                   
#' out <- GM(Y, list_hyper, list_init, track_q = T, track_vb_elbo = T)
#' @export
#'

GM  <- function(Y,
                list_hyper,
                tol = 1e-3,
                maxit = 1e5,
                #
                verbose = T,
                anneal = NULL,
                track_q = F,
                track_vb_elbo = F,
                set_v0 = 0.05) {
  # 
  pt <- Sys.time()
  
  
  list_init <- list(a_rho = 1, b_rho = 1, a_tau = 1, b_tau = 1)
                     
  
  args <-
    list(
      Y = Y,
      list_hyper = list_hyper,
      list_init = list_init,
      tol = tol,
      maxit = maxit,
      verbose = verbose,
      anneal = anneal,
      track_q = track_q,
      track_vb_elbo = track_vb_elbo
    )
  
  # dim
  #
  p <- ncol(Y)
  n <- nrow(Y)
  
  list_hyper$v0 = set_v0
  # Gathering initial parameters
  #
  Y <- scale(Y, center = TRUE, scale = FALSE) # center the data
  S <- crossprod(Y)
  list2env(list_init, envir = .GlobalEnv)
  list2env(list_hyper, envir = .GlobalEnv)
  
  
  if (!'a_rho' %in% names(list_init) | ('a_rho' %in% names(list_init) & is.null(list_init$a_rho))) {
    a_rho <- 1
  }
  
  if (!'b_rho' %in% names(list_init) | ('b_rho' %in% names(list_init) & is.null(list_init$b_rho))) {
    b_rho <- 1
  }
  
  if (!'a_tau' %in% names(list_init) | ('a_tau' %in% names(list_init) & is.null(list_init$a_tau))) {
    a_tau <- 1
  }
  
  if (!'b_tau' %in% names(list_init) | ('b_tau' %in% names(list_init) & is.null(list_init$b_tau))) {
    b_tau <- 1
  }
  
  
  #v0 <- get_m_v0(a_v0, b_v0)
  
  
  if (!'Omega' %in% names(list_init) | ('Omega' %in% names(list_init) & is.null(list_init$Omega))) {
    # if init Omega is not given, use the empirical
    Sigma <- (S + v0 * diag(p)) / n
    Omega <- as.matrix(Matrix::nearPD(solve(Sigma))$mat)
  }
  
  #
  
  m_tau <- 1 # get_m_tau(a_tau, b_tau) # for m_delta
  m_log_rho <- get_m_log_rho(a_rho, b_rho)
  m_log_one_minus_rho <- get_m_log_one_minus_rho(a_rho, b_rho)
  # % #
  
  # anneal
  #
  if (!is.null(anneal)) {
    vec_c <- get_annealing_ladder_(anneal, verbose = verbose)
  } else {
    vec_c <- NULL
  }
  
  # track objective
  #
  if (track_q) {
    vec_Q <- c()
    vec_Q_diff <- c()
  } else{
    vec_Q <- NA
    vec_Q_diff <- NA
  }
  
  # track objective in the VB
  #
  if (track_vb_elbo) {
    list_elbo <- list()
  } else{
    list_elbo <- NA
  }
  
  #
  eps <- .Machine$double.eps
  it <- 0
  n_warning <- 0
  vec_vb_it <- c()
  vec_n_warning_vb <- c()
  Q_diff <- Inf
  Q_old <- -Inf
  
  while (Q_diff > tol & (it < maxit)) {
    #
    it <- it + 1
    
    if (verbose != 0 & it %% 5 == 0 & T)
      cat(paste0("Iteration ", format(it), "... \n"))
    
    
    # VB step :
    # ====== #
    
    lb_diff <- Inf
    lb_old <- -Inf
    vb_it <- 0
    n_warning_vb <- 0
    
    if (track_vb_elbo) {
      tmp <- c()
    }
    
    while (lb_diff > tol & (vb_it < maxit)) {
      #
      vb_it <- vb_it + 1
      
      if (verbose != 0 & vb_it %% 5 == 0 & F)
        cat(paste0("VB iteration ", format(vb_it), "... \n"))
      
      if (!is.null(vec_c) && vb_it <= length(vec_c)) {
        vbc <- vec_c[vb_it]
        print(paste0("Temperature: ", format(1 / vbc, digits = 3)))
      } else {
        vbc <- 1
      }
      
      #v0 <- get_m_v0(a_v0, b_v0)
      # % #
      m_delta <- update_m_delta(Omega, m_tau, m_log_rho, m_log_one_minus_rho, v0, v1, vbc)
      E1 <- get_E1(m_delta, v0, v1)
      # % #
      # % # v0
      #a_v0 <- update_a_v0(av, m_delta)
      #b_v0 <- update_b_v0(bv, m_delta, Omega)
      #v0 <- get_m_v0(a_v0, b_v0)
      # % # tau
      a_tau <- update_a_tau(a, p, vbc)
      b_tau <- update_b_tau(Omega, E1, b, vbc)
      m_tau <-  1 # get_m_tau(a_tau, b_tau)
      m_log_tau <- 0 # get_m_log_tau(a_tau, b_tau)
      # % #
      
      # % # rho
      a_rho <- update_a_rho(m_delta, ar, vbc)
      b_rho <- update_a_rho(m_delta, br, vbc)
      m_log_rho <- get_m_log_rho(a_rho, b_rho)
      m_log_one_minus_rho <- get_m_log_one_minus_rho(a_rho, b_rho)
      # % #
      
      #
      vb_step_check <- T
      if (isTRUE(all.equal(vbc, 1))) {
        #
        lb <-  get_elbo_GM(Omega,
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
                           p
        )
        
        lb_diff <- abs(lb - lb_old)
        
        # check the increment in ELBO
        if (vb_step_check && lb + eps < lb_old) {
          warning(paste0(
            "Non-increasing in the ELBO: lb0 = ",
            lb_old,
            ", lb = ",
            lb,
            '\n'
          ))
          n_warning_vb <- n_warning_vb + 1
        }
        
        if (verbose != 0 & vb_it %% 5 == 0 & F)
          cat(paste0(
            "Difference lb from previous iteration: ",
            format(lb_diff),
            "\n"
          ))
        
        lb_old <- lb
        
        # track
        if (track_vb_elbo) {
          tmp <- c(tmp, lb)
        }
      }
    }
    if (track_vb_elbo) {
      list_elbo <- c(list_elbo, list(tmp))
    }
    
    #
    vec_vb_it <- c(vec_vb_it, vb_it)
    vec_n_warning_vb <- c(vec_n_warning_vb, n_warning_vb)
    
    if (vb_it == maxit) {
      warning('maxiter is reached. The VB algorithm may not converge ... \n')
    }
    
    # M step :
    # ====== #
    m_step_check <- TRUE
    
    # increment in ELBO in M-step
    
    if (m_step_check) {
      Q0 <- get_elbo_GM(Omega,
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
                        p
      )
    }
    
    # % # Omega
    bool_cpp <- F
    if (bool_cpp) {
      bool_direct_solve <- F # keep F as otherwise Omega inverted at each iteration
      if (bool_direct_solve) {
        out <- M_Omega_direct_solve(n, p, Omega, S, lambda, m_tau * E1)
        Omega <- out$Omega
      } else {
        out <- M_Omega(n, p, Sigma, Omega, S, lambda, m_tau * E1)
        Omega <- out$Omega
        Sigma <- out$Sigma
      }
    } else {
      Omega <- get_omega(m_tau * E1, S, Omega, lambda, n, p)
    }
    # % #
    
    Q <- get_elbo_GM(Omega,
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
                     p
    )
    
    Q_diff <- abs(Q - Q_old)
    
    if (m_step_check && Q + eps < Q0) {
      warning(paste0("Non-increasing in the M-step : Q0 = ", Q0, ", Q = ", Q))
      n_warning <- n_warning + 1
    }
    
    if (verbose != 0 & it %% 5 == 0 & F) ## note i have included F to turn this printing off
      cat(paste0(
        "Difference Q from previous iteration: ",
        format(Q_diff),
        "\n"
      ))
    
    Q_old <- Q
    
    if (track_q) {
      vec_Q <- c(vec_Q, Q)
      vec_Q_diff <- c(vec_Q_diff, Q_diff)
    }
  }
  
  if (it == maxit) {
    warning('maxiter is reached. The algorithm may not converge ... \n')
  }
  
  pt <- Sys.time() - pt
  #print(pt)
  
  create_named_list_(
    list_hyper,
    list_init,
    Omega,
    m_delta,
    E1,
    a_tau,
    b_tau,
    a_rho,
    b_rho,
    it,
    n_warning,
    vec_Q,
    vec_Q_diff,
    vec_vb_it,
    vec_n_warning_vb,
    list_elbo,
    pt,
    Y = Y,
    v0
  )
}
