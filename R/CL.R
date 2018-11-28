
is_acyclic <- function(adj_matrix){
  if (any(diag(adj_matrix) == 1)){
    stop("The graph represented by the matrix contains loops.")
  }
  while (any(rowSums(adj_matrix) == 1 |
             rowSums(adj_matrix) == 0)){
    idx <- which(rowSums(adj_matrix) == 1 |
                   rowSums(adj_matrix) == 0)[1]
    adj_matrix <- as.matrix(adj_matrix[- idx, - idx])
  }

  if (nrow(adj_matrix) == 0){
    res <- TRUE
  } else {
    res <- FALSE
  }
  res
}

ChowLiu <- function(data, root = NULL, ...){
  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors")
  }

  # Calculating mutual information
  nodes <- names(data)
  n_var <- length(nodes)

  pair_1 <- rep(nodes[- n_var], (n_var - 1):1)
  pair_2 <- unlist(sapply(1:(n_var - 1), function(n){
    nodes[- (1:n)]
    }))

  MI_fun <- function(var1, var2){
    MI2(data[, var1], data[, var2], ...)
  }

  MI <- mapply(MI_fun, pair_1, pair_2)
  MI_tab <- data.frame(pair_1, pair_2, MI)

  ord_idx <- order(MI_tab$MI, decreasing = TRUE)
  MI_tab <- MI_tab[ord_idx, ]
  rownames(MI_tab) <- NULL

  # Construct skeleton for Chow-Liu tree
  adj_matrix <- matrix(0, nrow = n_var, ncol = n_var)
  rownames(adj_matrix) <- colnames(adj_matrix) <- nodes
  i <- 1

  while (sum(adj_matrix) < 2 * (n_var - 1)){
    adj_matrix_temp <- adj_matrix
    edge_1 <- as.character(MI_tab[i, 1])
    edge_2 <- as.character(MI_tab[i, 2])
    adj_matrix_temp[edge_1, edge_2] <-
      adj_matrix_temp[edge_2, edge_1] <- 1

    if (is_acyclic(adj_matrix_temp)){
      adj_matrix <- adj_matrix_temp
    }

    i <- i+1
  }

  # Determine DAG
  if(is.null(root)){
    root <- sample(nodes,1)
  }

  adj_matrix_directed <- matrix(NA, nrow = n_var, ncol = n_var)
  rownames(adj_matrix_directed) <-
    colnames(adj_matrix_directed) <- nodes

  while (nrow(adj_matrix) > 0) {
    adj_matrix_directed[root, ] <- adj_matrix[root, ]
    adj_matrix <- adj_matrix[! nodes %in% root, ]
    adj_matrix[,root] <- 0
    #opdater rod
  }
}




