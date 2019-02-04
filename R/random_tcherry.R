
#' @export
random_tcherry <- function(n, n_levels, noise = NULL){
  if (length(n) != 1 | ! is.numeric(n)){
    stop("n must be a single integer.")
  }

  if (n %% 1 != 0 | n <= 1){
    stop("n must be a positive integer and at least 2.")
  }

  if (length(noise) != 1 | ! is.numeric(noise)){
    stop("noise must be a single numeric number.")
  }

  if (noise < 0){
    stop("noise must be non-negative.")
  }

  if (! is.vector(n_levels) | ! is.numeric(n_levels)){
    stop("n_levels must be a numeric vector.")
  }

  if (! all(n_levels %% 1 == 0) | ! all(n_levels >= 1)){
    stop("n_levels must be all positive integers.")
  }

  if (length(n_levels) != n){
    stop("There are not enough specified number of levels for the
         number of variables.")
  }



  var_names <- paste("V", 1:n, sep = "")
  adj_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(adj_matrix) <- colnames(adj_matrix) <- var_names

  idx_1_edge <- sample(1:n, 2)
  adj_matrix[idx_1_edge[1], idx_1_edge[2]] <-
    adj_matrix[idx_1_edge[2], idx_1_edge[1]] <- 1

  tcherry_nodes <- var_names[idx_1_edge]
  nodes_remaining <- setdiff(var_names, tcherry_nodes)

  level_names_1 <- paste("l", 1:n_levels[idx_1_edge[1]], sep = "")
  CPT1 <- as.table(stats::runif(n_levels[idx_1_edge[1]]))
  CPT1 <- CPT1 / sum(CPT1)
  dimnames(CPT1) <- list(level_names_1)
  names(dimnames(CPT1)) <- var_names[idx_1_edge[1]]

  level_names_2 <- paste("l", 1:n_levels[idx_1_edge[2]], sep = "")
  CPT2 <- as.table(stats::runif(n_levels[idx_1_edge[2]]))
  CPT2 <- CPT2 / sum(CPT2)
  dimnames(CPT2) <- list(level_names_2)
  names(dimnames(CPT2)) <- var_names[idx_1_edge[2]]

  CPTs <- list(CPT1, CPT2)
  names(CPTs) <- var_names[idx_1_edge]

  i <- 3
  while (length(nodes_remaining) != 0){
    new_var <- sample(nodes_remaining, 1)
    edges_idx <- which(adj_matrix == 1, arr.ind = TRUE)
    new_edge_idx <- edges_idx[sample(1:nrow(edges_idx), 1), ]
    edge_1_var <- var_names[new_edge_idx[1]]
    edge_2_var <- var_names[new_edge_idx[2]]

    adj_matrix[new_var, new_edge_idx[1]] <-
      adj_matrix[new_edge_idx[1], new_var] <-
      adj_matrix[new_var, new_edge_idx[2]] <-
      adj_matrix[new_edge_idx[2], new_var] <- 1
    nodes_remaining <- setdiff(nodes_remaining, new_var)

    l_new_var <- n_levels[var_names == new_var]
    l_edge_1 <- n_levels[new_edge_idx[1]]
    l_edge_2 <- n_levels[new_edge_idx[2]]
    CPTnew <- as.table(stats::runif(l_new_var * l_edge_1 * l_edge_2))
    dim(CPTnew) <- c(l_new_var, l_edge_1, l_edge_2)

    if (! is.null(noise)){
      for (j in 2:l_edge_2) {
        CPTnew[, , j] <- CPTnew[, , 1] + abs(rnorm(l_new_var * l_edge_1,
                                             sd = noise))
      }
    }

    cs <- colSums(CPTnew)
    CPTnew <- sweep(CPTnew, c(2, 3), cs, FUN = "/")
    levels1 <- paste("l", 1:l_new_var, sep = "")
    levels2 <- paste("l", 1:l_edge_1, sep = "")
    levels3 <- paste("l", 1:l_edge_2, sep = "")
    dimnames(CPTnew) <- list(levels1, levels2, levels3)
    names(dimnames(CPTnew)) <- c(new_var, edge_1_var, edge_2_var)
    CPTs[[i]] <- CPTnew
    names(CPTs)[i] <- new_var
    i <- i + 1
  }
  return(list("adj_matrix" = adj_matrix,
              "CPTs" = CPTs))
}
