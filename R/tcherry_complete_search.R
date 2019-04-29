#' Determine the k'th order t-cherry tree from data with the greatest
#' weight
#'
#' @description Determine the structure of the k'th order t-cherry tree
#' from data with the greatest weight based on a complete search.
#'
#' @param data The data the tree structure should be based on.
#' @param k The order of the t-cherry tree.
#' @param ... Additional arguments passed to \code{weight_junction_tree}.
#'
#' @details
#'
#' The algorithm for constructing the t-cherry tree from
#' data is based on an atempt to minimize the Kullback-Leibler
#' divergence, by mazimising the weight. All possible structures are
#' determined and the one with the highest weight is chosen.
#'
#' Note that this procedure is highly inefficient, and only suited for
#' small problems.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{model} A list with the following components:
#' \itemize{
#' \item \code{weight} The weight of the final k'th order t-cherry tree.
#' \item \code{cliques} A list containing the cliques of
#'  the k'th order t-cherry tree.
#' \item \code{separators} A list containing the separators of a
#' junction tree for the k'th order t-cherry tree.
#' \item \code{adj_matrix} The adjacency matrix for the k'th order
#' t-cherry tree.
#' \item \code{n_edges} The number of edges in the resulting graph.
#' }
#' \item \code{n_models} The number of considered models.
#' }
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @seealso \code{\link{weight_junction_tree}} for calculation of the
#' weight and \code{\link{k_tcherry_step}} for a more efficient but
#' greedy algorithm.
#'
#' @examples
#' set.seed(43)
#' var1 <- c(sample(c(1, 2), 100, replace = TRUE))
#' var2 <- var1 + c(sample(c(1, 2), 100, replace = TRUE))
#' var3 <- var1 + c(sample(c(0, 1), 100, replace = TRUE,
#'                         prob = c(0.9, 0.1)))
#' var4 <- c(sample(c(1, 2), 100, replace = TRUE))
#' var5 <- var2 + var3
#' var6 <- var1 - var4 + c(sample(c(1, 2), 100, replace = TRUE))
#' var7 <- c(sample(c(1, 2), 100, replace = TRUE))
#'
#' data <- data.frame("var1" = as.character(var1),
#'                    "var2" = as.character(var2),
#'                    "var3" = as.character(var3),
#'                    "var4" = as.character(var4),
#'                    "var5" = as.character(var5),
#'                    "var6" = as.character(var6),
#'                    "var7" = as.character(var7))
#'
#' (tch <- tcherry_complete_search(data, 5, smooth = 0.1))
#'
#' # For plotting
#' library(gRbase)
#' library(Rgraphviz)
#' tcherry_tree <- as(tch$model$adj_matrix, "graphNEL")
#' plot(tcherry_tree)
#'
#' # For probability propagation
#' library(gRain)
#' model <- grain(tcherry_tree, data = data, smooth = 0.1)
#' querygrain(model)
#' @export

tcherry_complete_search <- function(data, k, ...){

  if (any(is.na(data))){
    warning(paste("The data contains NA values.",
                  "Theese will be excluded from tables,",
                  "which may be problematic.",
                  "It is highly recommended to manually take",
                  "care of NA values before using the data as input.",
                  sep = " "))
  }

  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix.")
  }

  data <- as.data.frame(data)

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors.")
  }

  if (length(k) != 1){
    stop("k must be a single positive integer.")
  }

  if (k %% 1 != 0 | k <= 1){
    stop("k must be a positive integer and at least 2.")
  }

  nodes <- names(data)
  n_var <- length(nodes)

  first_cliques <- utils::combn(nodes, k)

  models <- apply(first_cliques, 2, function(col){
    adj_matrix <- matrix(0, nrow = n_var, ncol = n_var)
    rownames(adj_matrix) <- colnames(adj_matrix) <- nodes
    adj_matrix[col, col] <- 1
    diag(adj_matrix[col, col]) <- 0
    list("cliques" = list(col),
         "separators" = list(),
         "unused" = setdiff(nodes, col),
         "adj_matrix" = adj_matrix)
  })

  n_iter <- n_var - k

  if (n_iter != 0){
    for (iter in 1:n_iter) {
      new_models <- list()
      idx <- 1
      for (model in models) {
        for (clique in model$cliques) {
          seps <- utils::combn(clique, k - 1)
          seps <- split(seps, rep(1:ncol(seps), each = nrow(seps)))
          for (sep in seps) {
            for (var in model$unused) {
              new_model <- list()
              new_clique <- c(sep, var)

              new_matrix <- model$adj_matrix
              new_matrix[new_clique, new_clique] <- 1
              diag(new_matrix[new_clique, new_clique]) <- 0

              new_model$cliques <- c(model$cliques, list(new_clique))
              new_model$separators <- c(model$separators, list(sep))
              new_model$unused <- setdiff(model$unused, var)
              new_model$adj_matrix <- new_matrix

              new_models[[idx]] <- new_model
              idx <- idx + 1
            }
          }
        }
      }
      matrix_list <- lapply(new_models, function(l){
        l$adj_matrix
      })
      duplic <- duplicated(matrix_list)
      new_models <- new_models[! duplic]
      models <- new_models
    }
  }

  weights <- sapply(models, function(model){
    weight_junction_tree(model$cliques, model$separators, data, ...)
  })

  idx.max <- which.max(weights)
  model <- models[[idx.max]]
  model$weight <- weights[idx.max]
  model <- model[- 3]
  model$n_edges <- sum(model$adj_matrix) / 2

  return(list("model" = model,
              "n_models" = length(models)))

}

