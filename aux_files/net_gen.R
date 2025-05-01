


### create_c_networks returns Cn adjacency matrices (networks) and Cn corresponding datasets. 
### The adjacency matrices have dimension PxP and are all scale-free.
### Sparsity of the adjacency matrices are controlled by the 'sparsity' argument.
### Seed can be set in the argument of the function

### outputs: $Omegas (list of simulated precision matrices), $Ys (list of generated datasets),
### $Ns (list of number of rows in each Y of Ys), $pos_id (indices of increasing entries in the
### precision matrices), $neg_id (indices of increasing entries in the precision matrices), 
### $As (Adjacency matrices corresponding to $Omegas)

create_Cn_networks <- function(P = 50, Cn = 3, seed_1 = 123, seed_2 = 123, 
                              frac_change = 0.2, Ns_sample = c(150,150,150)){
  
  
  
  set.seed(seed_1)
  
  # Sample sizes for each dataset
  Ns <- sample(Ns_sample, Cn, replace = T)
  
  #generate an initial network
  net0 <- huge::huge.generator(n = Ns[1], d = P, graph = 'scale-free', v = 0.4, u = 0.05, verbose = F)
  diag.omega <- diag(net0$omega)
  net0$omega[as.numeric(net0$theta) != 1] <- 0
  diag(net0$omega) <- diag.omega
  
  pol <- matrix(data = sample(c(-1,1), P^2, replace = T), nrow = P)
  pol <- upper.tri(pol)*pol + t(upper.tri(pol))*t(pol)
  diag(pol) <- 1
  net0$omega <- net0$omega * pol
  
  # Replicating the initial precision matrix. 
  Omegas <- replicate(Cn, net0$omega, simplify=FALSE)
  
  # Assert which vertices have an edge (edge_id) between them and which don't (non_edge_id)
  edge_id <- which(as.logical(net0$theta)  & upper.tri(net0$omega), arr.ind = T)
  non_edge_id <- which(!as.logical(net0$theta) & upper.tri(net0$omega), arr.ind = T)
  
  
  npos <- nneg <- ceiling(nrow(edge_id)*frac_change)
  
  
  # Of the absent edges, randomly select 'npos' of them to have a positive correlation with the covariate
  # Of the absent edges, randomly select 'nneg' of them to have a negative correlation with the covariate
  neg_id <- matrix(edge_id[sample(1:nrow(edge_id),nneg),], ncol=2)
  pos_id <- matrix(non_edge_id[sample(1:nrow(non_edge_id),npos),], ncol = 2)
  


    
  # Iterate over each disappearing edge and have them 'fade' out linearly
  for (i in 1:nneg) {
    val <- net0$omega[neg_id[i,1],neg_id[i,2]]
    #k <- sample(2:Cn, 1, replace = T)
    for (t in 2:Cn) {
    Omegas[[t]][neg_id[i,1],neg_id[i,2]] <- Omegas[[t]][neg_id[i,2],neg_id[i,1]] <- ((Cn-t)/(Cn-1))*val
    }
  }
    
  # Iterate over each appearing edge and have them 'fade' in linearly
  hi_val <- max(abs(net0$omega[lower.tri(net0$omega)]))
  for (i in 1:npos) {
    #k <- sample(2:Cn, 1, replace = T)
    val <- hi_val - sample(c(0,0.02,0.04,0.06,0.08,0.1), 1, replace = T)
    plus_min <- sample(c(-1,1),1,replace = T)
    for (t in 2:Cn) {
      Omegas[[t]][pos_id[i,1],pos_id[i,2]] <- Omegas[[t]][pos_id[i,2],pos_id[i,1]] <- plus_min*val*(t-1)/(Cn-1)
    }
  }
  
  
  Ys = list()
  
  for(i in seq_along(Omegas)){
    Omegas[[i]] <- (Omegas[[i]] + t(Omegas[[i]]))/2
  }
  for(i in seq_along(Omegas)){
    while(!matrixcalc::is.positive.definite(Omegas[[i]])){
     nets <- create_Cn_networks(P = P, Cn = Cn, seed_1 = seed_1 + 1, seed_2 = 123,
                                frac_change = frac_change)
     return(nets)
    }
  }
  set.seed(seed_2)

  for(y in 1:Cn){
    
    Ys[[y]] <- MASS::mvrnorm(n=Ns[y], mu = rep(0,P), Sigma = solve(Omegas[[y]]))
    Ys[[y]] <- scale(Ys[[y]])
    
  }
  #extract adjacencies
  
  As <- lapply(Omegas, function(omega) (abs(omega) > 0) * 1)
  
  cat("seed that works for P =",P,", frac =",frac_change," : ",seed_1,"\n")

  return(list(Omegas = Omegas, Ys = Ys, Ns = Ns, neg_id = neg_id, As = As,
              net0_sig = net0$sigma, pos_id = pos_id, seed_1 = seed_1, edge_id = edge_id))

}


plot_adjacency_matrices <- function(adj_list) {
  # Check if input is a list
  if (!is.list(adj_list)) stop("Input must be a list of adjacency matrices.")
  
  # Iterate over the adjacency matrices
  for (i in seq_along(adj_list)) {
    # Access the matrix and ensure it's numeric
    adj_matrix <- as.matrix(adj_list[[i]])
    if (!is.matrix(adj_matrix) || !is.numeric(adj_matrix)) {
      warning(paste("Element", i, "is not a valid numeric matrix. Skipping..."))
      next
    }
    
    # Set the diagonal to 0
    diag(adj_matrix) <- 0
    
    # Convert to an igraph object
    graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
    
    # Plot the graph
    plot(graph, main = paste("Graph", i))
  }
}


