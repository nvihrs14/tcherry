#' Test whether an undirected graph is acyclic
#'
#' @description Test whether an undirected graph, represented by
#' an adjacency matrix, is acyclic.
#'
#' @param adj_matrix The adjacency matrix representing the graph.
#' @details Notice that the function cannot cope with loops.
#' If the graph has loops, an error is returned.
#' @return A logical value indicating whether the graph is acyclic.
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#' @examples
#' adj_matrix_cyclic <- matrix(c(0, 1, 1, 1,
#'                               1, 0, 0, 1,
#'                               1, 0, 0, 0,
#'                               1, 1, 0, 0),
#'                              nrow = 4)
#'
#' is_acyclic(adj_matrix_cyclic)
#'
#' adj_matrix_acyclic <- matrix(c(0, 0, 1, 1,
#'                                0, 0, 0, 1,
#'                                1, 0, 0, 0,
#'                                1, 1, 0, 0),
#'                               nrow = 4)
#'
#' is_acyclic(adj_matrix_acyclic)
#' @export

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

#' Determines the Chow-Liu tree for data
#'
#' @description Determines the structure and the conditional
#' probability tables for the Chow-Liu tree fitted to data.
#'
#' @param data The data set the model should be fitted for.
#' @param root An optional argument, choosing a userspecified
#' root for the tree.
#' @param bayes_smooth Additional cell counts for bayesian
#' estimation.
#' @param ... Additional parameters passed to \code{MI2}.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{skeleton_adj} The adjacency matrix for the skeleton
#' of the Chow-Liu tree.
#' \item \code{adj_DAG} The adjacency matrix of the resulting DAG.
#' \item \code{CPTs} The estimated conditional probability tables
#' of the bayesian network.
#' }
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @seealso \code{\link{MI2}} for mutual information,
#' \code{\link{CPT}} for probability tables and
#' \code{\link{is_acyclic}} for a test for cycles.
#'
#' @examples
#' set.seed(43)
#' var1 <- c(sample(c(1, 2), 50, replace = TRUE))
#' var2 <- var1 + c(sample(c(1, 2), 50, replace = TRUE))
#' var3 <- var1 + c(sample(c(0, 1), 50, replace = TRUE,
#'                         prob = c(0.9, 0.1)))
#' var4 <- c(sample(c(1, 2), 50, replace = TRUE))
#'
#' data <- data.frame("var1" = as.character(var1),
#'                    "var2" = as.character(var2),
#'                    "var3" = as.character(var3),
#'                    "var4" = as.character(var4))
#'
#' CL <- ChowLiu(data, root = 'var1', smooth = 0.1)
#' @export

ChowLiu <- function(data, root = NULL, bayes_smooth = 0, ...){
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

  if (! root %in% nodes){
    stop("The specified root is not a node.")
  }

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

  skeleton_adj <- adj_matrix

  # Determine DAG
  if(is.null(root)){
    root <- sample(nodes,1)
  }

  adj_matrix_directed <- matrix(NA, nrow = n_var, ncol = n_var)
  rownames(adj_matrix_directed) <-
    colnames(adj_matrix_directed) <- nodes

  while (nrow(adj_matrix) > 0) {
    adj_matrix_directed[root, ] <- adj_matrix[root, ]
    adj_matrix <- adj_matrix[! rownames(adj_matrix) %in% root, ,
                             drop = FALSE]
    adj_matrix[, root] <- 0
    kids_idx <- which(colSums(adj_matrix_directed[root, ,
                                      drop = FALSE]) > 0)
    root <- nodes[kids_idx]
  }

  # Calculate conditional probability tables
  CPTs <- CPT(adj_matrix_directed, data,
              bayes_smooth = bayes_smooth)

  return(list("skeleton_adj" = skeleton_adj,
              "adj_DAG" = adj_matrix_directed,
              "CPTs" = CPTs))
}




