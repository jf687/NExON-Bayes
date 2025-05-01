## 'evaluate_network' takes a single true and estimated network and returns the confusion matrix components

evaluate_network <- function(true_network, estimated_network,t) {
  if (!all(dim(true_network) == dim(estimated_network))) {
    stop("The dimensions of the true network and the estimated network must match.")
  }
  
  
  
  LT <- lower.tri(true_network)
  TP <- sum(true_network[LT] >t & estimated_network[LT] > t)
  TN <- sum(true_network[LT] <t & estimated_network[LT] < t)
  FP <- sum(true_network[LT] <t & estimated_network[LT] > t)
  FN <- sum(true_network[LT] >t & estimated_network[LT] < t)
  
  return(list(
    true_positives = TP,
    true_negatives = TN,
    false_positives = FP,
    false_negatives = FN,
    conf.mat = matrix(c(TP, FN, FP, TN), nrow = 2, byrow = T)
  ))
}

## 'evaluate_network_list' takes multiple true and estimated networks and returns the confusion matrix components
## It also returns the confusion matrix in itself, and 'P' which is the dimension of the networks that it is calculating.

evaluate_network_list <- function(true_networks, estimated_networks,t) {
  # Ensure both lists have the same length
  if (length(true_networks) != length(estimated_networks)) {
    stop("The lists of true networks and estimated networks must have the same length.")
  }
  
  total_TP <- 0
  total_TN <- 0
  total_FP <- 0
  total_FN <- 0
  
  conf.mat <- matrix(0, nrow = 2, ncol = 2)
  
  P <- ncol(true_networks[[1]])
  

  for (i in seq_along(true_networks)) {
    # Call the original evaluate_network function
    result <- evaluate_network(true_networks[[i]], estimated_networks[[i]],t )

    total_TP <- total_TP + result$true_positives
    total_TN <- total_TN + result$true_negatives
    total_FP <- total_FP + result$false_positives
    total_FN <- total_FN + result$false_negatives
    conf.mat <- conf.mat + result$conf.mat
  }
  

  return(list(
    P = P,
    true_positives = total_TP,
    true_negatives = total_TN,
    false_positives = total_FP,
    false_negatives = total_FN,
    conf.mat = conf.mat)
  )
}

## 'precision' gives the precision of estimations based on a confusion matrix. 
## precision = TN/(TN + )

recall = function (conf.mat) {

  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  }
  else if (conf.mat[1, 1] == 0 & conf.mat[1, 2] == 0) {
    return(1)
  }
  else {
    return(conf.mat[1, 1]/(conf.mat[1, 1] + conf.mat[1, 2]))
  }
}

precision = function (conf.mat) {

  if (conf.mat[1, 1] == 0 & conf.mat[2, 1] == 0) {
    return(1)
  }
  else {
    return(conf.mat[1, 1]/(conf.mat[1, 1] + conf.mat[2, 1]))
  }
}

sparsity <- function(omega, strict = TRUE, threshold = 1e-3){
  # Calculate the total number of elements in the matrix
  
  diag(omega) = 0
  p=ncol(omega)
  
  total_elements <- length(omega)
  
  # Calculate the number of zero elements in the matrix, or the number below a given threshold
  if (strict){
    zero_elements <- sum(omega == 0)
  } else {
    zero_elements <- sum(omega < threshold)
  }
  
  # Calculate the sparsity as 1 minus the proportion of zero elements
  #sparsity <- 1 - zero_elements / total_elements
  sparsity =  1 - (zero_elements -p)/ (total_elements-p)
  
  return(sparsity)
}

MCC = function(g,g.hat){
  p = nrow(g[,])
  diag(g) = rep(0,p) # Let diagonal elements be zero
  diag(g.hat) = rep(0,p) 
  tp = sum(g.hat ==1 & g ==1)/10 # True positives. Divide by 10 to avoid integer overflow. 
  fp = sum(g.hat ==1 & g ==0)/10 # False positives
  tn = (sum(g.hat == 0 & g == 0) - p)/10 # True negatives (do not include diagonal elements)
  fn = sum(g.hat == 0 & g == 1)/10 # False negatives
  return((tp*tn - fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))))
}

plot_confusion_matrix <- function(conf.mat){
  True <- factor(c(1, 1, 0, 0))
  Predicted <- factor(c(1, 0, 1, 0))
  conf.mat <- as.vector(t(conf.mat))
  df <- data.frame(True, Predicted, conf.mat)

  ggplot(data =  df, mapping = aes(x = True, y = Predicted)) +
    geom_tile(aes(fill = conf.mat), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", conf.mat)), vjust = 1) +
    scale_fill_gradient(low = "white", high = "lightblue") +
    theme_bw() + theme(legend.position = "none")
}

DKL <- function(Om.est, Om.true){
  # if(dim(Om.est) == dim(Om.true)){
  #   stop("Discrepant matrix dimensions.")
  # }
  om.om <- solve(Om.true) %*% Om.est
  det.om.om <- determinant(om.om, logarithm = T)$modulus[[1]]
  P <- ncol(Om.est)
  return( sum(diag(om.om))- det.om.om - P )
}

sparsity <- function(network){
  return(sum(network[lower.tri(network)] != 0)/sum(lower.tri(network)))
}