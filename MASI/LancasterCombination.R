spearman.matrix <- function(p.values){
  methods.names <- colnames(p.values)
  spearman.correlation <- matrix(0, nrow = ncol(p.values), ncol = ncol(p.values))
  rownames(spearman.correlation) <- methods.names
  colnames(spearman.correlation) <- methods.names
  for (i in 1:ncol(spearman.correlation)) {
    temp.i <- p.values[,i]
    for (j in 1:ncol(spearman.correlation)) {
      temp.j <- p.values[,j]
      spearman.correlation[i,j] <- cor(x = temp.i, y = temp.j, method = c("spearman"))
    }
  }
  diag(spearman.correlation) <- 0
  spearman.correlation[spearman.correlation < 0] <- 0
  return(spearman.correlation)
}

lancaster.combination <- function(Pvals, weight = TRUE, trimmed = 0.2){
  
  #### check if there is NA
  if (sum(is.na(Pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals < 0)+sum(Pvals > 1)) > 0){
    stop("P-values must be between 0 and 1!")
  }
  
  
  ####### check the validity of the user supplied weights and standadize them.
  if(weight){
    weight.matrix <- spearman.matrix(p.values = Pvals)
    temp <- Matrix::rowSums(weight.matrix)/ncol(weight.matrix)
    weight.pval <- temp
  }else{
    num_pval <- ncol(Pvals)
    weight.pval <- rep(2, num_pval)
    names(weight.pval) <- colnames(Pvals)
  }
  ###### calculate combanation p.values for each gene
  combination.pvalues.lancaster <- function(x, weight){
    location <- x == 1
    temp <- sort(x[! location], decreasing = FALSE)
    num <- floor(length(temp) * trimmed)
    Pvals.trimmed <- temp[(num + 1) : (length(temp) - num)]
    weight.use <- weight[names(Pvals.trimmed)]
    weight.use.norm <- weight.use/sum(weight.use)
    ####  Replace extreme p-values with 10e-320 to obtain an upper bound for the aggregated p-value. ####
    if (any(Pvals.trimmed < 9.99988867182683e-320)) {
      Pvals.trimmed[Pvals.trimmed < 9.99988867182683e-320] <- 9.99988867182683e-320
    }
    pval <- aggregation::lancaster(pvalues = Pvals.trimmed, weights = weight.use.norm)
  }
  pvals <- apply(Pvals, 1, combination.pvalues.lancaster, weight = weight.pval)
  combination.Pvals <- pvals
  combination.Pvals <- combination.Pvals[order(combination.Pvals)]
  return(data.frame(combination.Pvals))
}
