#' Gaussian graphical spike-and-slab models with time. 
#'
#' This function conducts simultaneous inference of a Gaussian graphical spike-and-slab models over time. This function allows varying intercepts across edges.
#'
#' @param Ys A list of data matrix of dimensions Nt x P, where Nt is the number of samples at time t and P is the number of nodes in the graph.
#' @param ts A vector of time points. 
#' 
#' @examples
#' Ys <- list()
#' set.seed(1234)
#' Nt1 <- 200
#' Nt2 <- 100
#' P <- 50
#' Ys[[1]] <- matrix(rnorm(Nt1 * P), Nt1, P)
#' Ys[[2]] <- matrix(rnorm(Nt2 * P), Nt2, P)
#' ts <- 1:2
#' ans <- NExON(Ys, ts, debug = T)

source('NExON/functions.R')

NExON <- function(Ys,
                ts,
                list_hyper = NULL,
                list_init = NULL,
                tol = 1e-3,
                maxit = 1e5,
                verbose = T,
                debug = F,
                v0_list = NULL
) {
  
  
  # Save inputs
  #
  args <-
    list(
      Ys = Ys,
      ts = ts,
      list_hyper = list_hyper,
      list_init = list_init,
      tol = tol,
      maxit = maxit,
      verbose = verbose,
      debug = debug
    )
  
  # Time
  #
  pt <- Sys.time()
  
  # Dimension
  #
  Ns <- sapply(Ys, nrow)
  
  if (length(unique(sapply(Ys, ncol))) != 1){
    stop("Same number of nodes should be considered across time points.")
  }else{
    P <- unique(sapply(Ys, ncol))
  }
  Tn <- length(Ns)
  
  
  #
  Ss <- lapply(Ys,crossprod)
  
  #
  if (verbose) cat("== Preparing the hyperparameters ... \n\n")
  
  list_hyper <- set_default(list_hyper, 'lambda', 2)
  if(list_hyper$lambda <= 0)stop("lambda must be positive.")
  
  if(is.null(v0_list)){
    v0_list <- rep(0.02, Tn)
  }
  
  for(v0 in v0_list){
    if(v0 <= 0)stop("All v0s must be positive.")
  }
  
  list_hyper <- set_default(list_hyper, 'v1', 100)
  if(list_hyper$v1 <= 0)stop("v1 must be positive.")
  if(any(v0_list >= list_hyper$v1))stop("v1 should be much greater than v0.")
  
  
  list_hyper <- set_default(list_hyper, 'a_sigma', 0.1)
  if(list_hyper$a_sigma <= 0)stop("a_sigma must be positive.")
  
  list_hyper <- set_default(list_hyper, 'b_sigma', 2)
  if(list_hyper$b_sigma <= 0)stop("b_sigma must be positive.")
  
  list_hyper <- set_default(list_hyper, 'n0', -2)
  
  list_hyper <- set_default(list_hyper, 't02', 5)
  if(list_hyper$t02 <= 0)stop("t02 must be positive.")
  
  if (verbose) cat("... done. == \n\n")
  
  if (verbose) cat("== Preparing the parameter initialisation ... \n\n")
  
  # Specify initial parameters unless provided
  # 
  lambda <- list_hyper$lambda
  v1 <- list_hyper$v1
  a_sigma <- list_hyper$a_sigma
  b_sigma <- list_hyper$b_sigma
  n0 <- list_hyper$n0
  t02 <- list_hyper$t02
  
  v0s <- v0_list
  #
  if(length(v0_list) != Tn){stop(paste0(length(v0_list)," values for v0 provided. Should be ",Tn))}
  #
  if (!'Omegas' %in% names(list_init) | ('Omegas' %in% names(list_init) & is.null(list_init$Omegas))) {
    
    Sigmas <- lapply(1:Tn, function(t)(Ss[[t]] + v0s[[t]] * diag(P))/Ns[t])
    list_init$Omegas <- lapply(Sigmas, function(x)as.matrix(Matrix::nearPD(solve(x))$mat))
    
  }else{
    
    for(t in 1:Tn) {
      if(nrow(list_init$Omegas[[t]])!=ncol(list_init$Omegas[[t]])){
        stop(paste0("Omega should be initialised as a square matrix at time ", t))
      }else if(!matrixcalc::is.symmetric.matrix(list_init$Omegas[[t]])){
        stop(paste0("Omega should be initialised as a symmetric matrix at time ", t))
      }else if(!matrixcalc::is.positive.definite(list_init$Omegas[[t]])){
        stop(paste0("Omega should be initialised as a positive definite matrix at time ", t))
      }
    }
    
  }
  
  if (!'mu_beta' %in% names(list_init) | ('mu_beta' %in% names(list_init) & is.null(list_init$mu_beta))) {
    
    list_init$mu_beta <- matrix(0, P, P)
    
  }else{
    
    # will only use the upper triangle 
    if(nrow(list_init$mu_beta) != ncol(list_init$mu_beta)){
      stop("mu_beta is not initialised as a square matrix.")
    }else if(nrow(list_init$mu_beta) != P){
      stop("Dimension of mu_beta does not match the dimension of precision matrices.")
    }
    
  }
  
  if (!'sig2_inv_beta' %in% names(list_init) | ('sig2_inv_beta' %in% names(list_init) & is.null(list_init$sig2_inv_beta))) {
    
    list_init$sig2_inv_beta <- matrix(1, P, P)
    
  }else{
    
    if(nrow(list_init$sig2_inv_beta) != ncol(list_init$sig2_inv_beta)){
      stop("sig2_inv_beta is not initialised as a square matrix.")
    }else if(nrow(list_init$sig2_inv_beta) != P){
      stop("Dimension of sig2_inv_beta does not match the dimension of precision matrices.")
    }else if(any(list_init$sig2_inv_beta <= 0 )){
      stop("sig2_inv_beta must be positive.")
    }
    
  }
  
  if (!'mu_zeta' %in% names(list_init) | ('mu_zeta' %in% names(list_init) & is.null(list_init$mu_zeta))) {
    
    list_init$mu_zeta <- matrix(list_hyper$n0, P, P)
    
  }else{
    
    # will only use the upper triangle 
    if(nrow(list_init$mu_zeta) != ncol(list_init$mu_zeta)){
      stop("mu_zeta is not initialised as a square matrix.")
    }else if(nrow(list_init$mu_zeta) != P){
      stop("Dimension of mu_zeta does not match the dimension of precision matrices.")
    }
    
  }
  
  if (!'sig2_inv_zeta' %in% names(list_init) | ('sig2_inv_zeta' %in% names(list_init) & is.null(list_init$sig2_inv_zeta))) {
    
    list_init$sig2_inv_zeta <-  matrix(1/list_hyper$t02,  P, P)
    
  }else{
    
    if(nrow(list_init$sig2_inv_zeta) != ncol(list_init$sig2_inv_zeta)){
      stop("sig2_inv_zeta is not initialised as a square matrix.")
    }else if(nrow(list_init$sig2_inv_zeta) != P){
      stop("Dimension of sig2_inv_zeta does not match the dimension of precision matrices.")
    }else if(any(list_init$sig2_inv_zeta <= 0 )){
      stop("sig2_inv_zeta must be positive.")
    }
    
  }
  
  
  if (!'alpha_sigma' %in% names(list_init) | ('alpha_sigma' %in% names(list_init) & is.null(list_init$alpha_sigma))) {
    
    list_init$alpha_sigma <- 1
    
  }else{
    
    if( list_init$alpha_sigma <= 0){
      stop("alpha_sigma must be positive.")
    }
    
  }
  
  if (!'beta_sigma' %in% names(list_init) | ('beta_sigma' %in% names(list_init) & is.null(list_init$beta_sigma))) {
    
    list_init$beta_sigma <- 1
    
  }else{
    
    if( list_init$beta_sigma <= 0){
      stop("beta_sigma must be positive.")
    }
    
  }

  
  
  
  #
  Omegas <- list_init$Omegas
  mu_zeta <- list_init$mu_zeta
  sig2_inv_zeta <- list_init$sig2_inv_zeta
  mu_beta <- list_init$mu_beta
  sig2_inv_beta <- list_init$sig2_inv_beta
  alpha_sigma <- list_init$alpha_sigma
  beta_sigma <- list_init$beta_sigma
  #
  
  # Initialise deduced quantities
  #

  m1_alphas <- get_m1_alphas(mu_zeta, mu_beta, ts)
  m2_beta <- get_m2_beta(mu_beta, sig2_inv_beta)
  
  if (verbose) cat("... done. == \n\n")
  
  # debug mode
  #
  if (debug) {
    
    # track ELBO within each variational step
    list_ELBO <- list()
    # track ELBO after each maximisation step
    vec_ELBO_CM <- c()
    
    # record number of warnings
    n_warning <- 0
    vec_n_warning_VB <- c()
    
  }
  
  #
  eps <- .Machine$double.eps
  
  #
  it <- 0
  vec_VB_it <- c()
  ELBO_CM_diff <- Inf
  ELBO_CM_old <- -Inf
  
  #
  while ((ELBO_CM_diff > tol & (it < maxit))) {
    
    # Iteration
    #
    it <- it + 1
    
    #if (verbose & it %% 5 == 0)
      #cat(paste0("Iteration ", format(it), "... \n"))
    
    
    # VBE step :
    # ======== #
    
    ELBO_diff <- Inf
    ELBO_old <- -Inf
    VB_it <- 0
    
    if(debug){
      n_warning_VB <- 0
      tmp <- c()
    }
    
    while (ELBO_diff > tol & VB_it < maxit) {
      
      #
      VB_it <- VB_it + 1
      
      #if (verbose != 0 & VB_it %% 5 == 0)
        #cat(paste0("VBE iteration ", format(VB_it), "... \n"))
      
      # % # m_deltas
      m_deltas <- update_m_deltas(Omegas, m1_alphas, v0s, v1, Tn) 
      dstars <- get_dstars(m_deltas, v0s, v1)
      # % #
      
      # % # z
      m_zs <- get_m_zs(m_deltas, m1_alphas, Tn)
      m2_zs <- get_m2_zs(m_zs, m1_alphas, Tn)
      
      # % # sigma
      alpha_sigma <- update_alpha_sigma(a_sigma, P)
      beta_sigma <- update_beta_sigma(m2_beta, b_sigma)
      
      m_sig2_inv <- get_m_sig2_inv(alpha_sigma, beta_sigma)
      m_log_sig2_inv <- get_m_log_sig2_inv(alpha_sigma, beta_sigma)
      # % #
      
      # % # zeta
      sig2_inv_zeta <- update_sig2_inv_zetaij(t02, P, Tn)
      mu_zeta <-
        update_mu_zetaij(mu_beta, sig2_inv_zeta, m_zs,ts, n0, t02, P)
      m2_zeta <- get_m2_zeta(mu_zeta, sig2_inv_zeta)
      # % #
      
      # % # beta
      sig2_inv_beta <- update_sig2_inv_beta(m_sig2_inv, ts, P)
      mu_beta <- update_mu_beta(sig2_inv_beta, m_zs, mu_zeta, ts, Tn)
      m2_beta <- get_m2_beta(mu_beta, sig2_inv_beta)
      # % #
      
      # % #
      m1_alphas <- get_m1_alphas(mu_zeta, mu_beta, ts)
      var_alphas <- get_var_alphas(sig2_inv_zeta, sig2_inv_beta, ts)
      # % #
      
      # % #
      ELBO <- get_elbo_zetaij(Omegas,
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
                              Tn)
      
      ELBO_diff <- abs(ELBO - ELBO_old)
      
      # check the increasing ELBO in the variational step
      #
      if (debug && ELBO + eps < ELBO_old) {
        
        warning(paste0(
          "Non-increasing in the VB step: ELBO_old = ",
          ELBO_old,
          ", ELBO = ",
          ELBO,
          '\n'
        ))
        n_warning_VB <- n_warning_VB + 1
      }
      
      if (verbose & VB_it %% 5 == 0)
        #cat(paste0(
        #  "Difference ELBO from previous iteration: ",
        #  format(ELBO_diff),
        #  "\n"
        #))
      
      ELBO_old <- ELBO
      
      if (debug) {
        tmp <- c(tmp, ELBO)
      }
      # % #
    }
    
    
    # % #
    if (debug) {
      list_ELBO <- c(list_ELBO, list(tmp))
      vec_n_warning_VB <- c(vec_n_warning_VB, n_warning_VB)
    }
    
    vec_VB_it <- c(vec_VB_it, VB_it)
    
    if (VB_it == maxit) {
      warning('Maximal number of iterations reached before convergence in VBE step. \n')
    }
    # % #
    
    # CM step :
    # ====== #
    #  check the increasing ELBO in the CM-step
    #
    if (debug) {
      ELBO_CM0 <- get_elbo_zetaij(Omegas,
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
                                  Tn)
      
    }
    

    Omegas <- get_omegas(lapply(dstars,function(x) x), 
                         Ss, Omegas, lambda, Ns, P, Tn)
    
    
    ELBO_CM <- get_elbo_zetaij(Omegas,
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
                               Tn)
    
    ELBO_CM_diff <- abs(ELBO_CM -   ELBO_CM_old)
    
    if (debug &&   ELBO_CM + eps <   ELBO_CM0) {
      
      warning(paste0("Non-increasing in the CM-step :   ELBO_0 = ",   ELBO_CM0, ",   ELBO = ",   ELBO_CM))
      n_warning <- n_warning + 1
      
    }
    
    if (verbose & it %% 5 == 0)
      cat(paste0(
        "Difference ELBO_CM from previous iteration: ",
        format(ELBO_CM_diff),
        "\n"
      ))
    
    ELBO_CM_old <- ELBO_CM
    
    if (debug) {
      vec_ELBO_CM <- c(vec_ELBO_CM, ELBO_CM)
    }
  }
  
  
  if(ELBO_CM_diff <= tol){
    if(verbose)
      cat(paste0("Convergence obtained after ", format(it), " iterations. \n",
                 "Optimal marginal log-likelihood variational lower bound ",
                 "(ELBO) = ", format(ELBO_CM), ". \n\n"))
  }
  
  if (it == maxit) {
    warning('Maximal number of iterations reached before convergence. Exit.')
  }
  
  pt <- Sys.time() - pt
  cat('Algorithm runtime: ',format(pt), '\n')
  
  estimates <- list( Omegas = Omegas,
                     m_deltas = m_deltas,
                     mu_beta = mu_beta,
                     sig2_inv_beta = sig2_inv_beta,
                     mu_zeta = mu_zeta,
                     sig2_inv_zeta = sig2_inv_zeta,
                     alpha_sigma = alpha_sigma,
                     beta_sigma = beta_sigma,
                     Ss = Ss # for model comparison
  )
  
  if(debug){
    
    debugs <- list( n_warning = n_warning,
                    vec_n_warning_VB = vec_n_warning_VB,
                    list_ELBO = list_ELBO,
                    vec_ELBO_CM = vec_ELBO_CM)
    
  }else{
    
    debugs <- NA
    
  }
  
  create_named_list_(
    args,
    estimates,
    debugs,
    it,
    vec_VB_it,
    pt,
    Ns,
    Ys = Ys
  )
}
