#' Determine the (k + 1)'th order t-cherry tree from a k'th order
#' t-cherry tree with the greatest weight
#'
#' @description Determine the structure of the (k + 1)'th order t-cherry
#' tree from a k'th order t-cherry tree with the greatest weight based
#' on a complete search.
#'
#' @param tch_cliq A list containing the cliques of a k'th order
#' t-cherry tree.
#' @param data The data the structure of the tree should be based on.
#' @param ... Additional arguments passed to \code{weight_junction_tree}.
#'
#' @details
#'
#' The algorithm for constructing the (k + 1)'th order t-cherry tree from
#' a k'th order t-cherry tree is based on an atempt to minimize the
#' Kullback-Leibler divergence, by mazimising the weight. All possible
#' structures are determined and the one with the highest weight is
#' chosen.
#'
#' Note that this procedure is highly inefficient, and only suited for
#' small problems.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{model} A list containing the following components:
#' \itemize{
#' \item \code{weight} The weight of the final (k + 1)'th order
#' t-cherry tree.
#' \item \code{cliques} A list containing the cliques of
#'  the (k + 1)'th order t-cherry tree.
#' \item \code{separators} A list containing the separators of a
#' junction tree for the (k + 1)'th order t-cherry tree.
#' \item \code{adj_matrix} The adjacency matrix for the (k + 1)'th order
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
#' weight and \code{\link{increase_order2}} for a more
#' efficient, but greedy algorithm.
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
#' ChowLiu_cliques <- list(c("var1", "var5"),
#'                         c("var2", "var5"),
#'                         c("var3", "var5"),
#'                         c("var3", "var7"),
#'                         c("var4", "var6"),
#'                         c("var5", "var6"))
#'
#' (tch <- increase_order_complete_search(ChowLiu_cliques, data,
#'                                        smooth = 0.1))
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

increase_order_complete_search <- function(tch_cliq, data, ...){

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

  if (! is.list(tch_cliq)){
    stop(paste("Cliques must be given in a list, each entry containing",
             "a vector with the names of the variables in the clique.",
             collapse = " "))
  }

  if (! compare::compare(unique(unlist(tch_cliq)), colnames(data),
                         ignoreOrder = TRUE)$result){
    stop(paste("The column names of data must be the same as the",
             "variable names in tch_cliq. All variables in data must",
             "be in at least one clique.", collapse = " "))
  }

  if (length(unique(sapply(tch_cliq, length))) != 1){
    stop(paste("tch_cliq should be the cliques of a k'th order t-cherry",
             "tree. Therefore they should all have the same length k.",
             collapse = " "))
  }

  # Reconstruct k and k'th order t-cherry tree.

  nodes <- names(data)
  n_var <- length(nodes)

  n_cliq <- length(tch_cliq)
  k <- length(tch_cliq[[1]])

  if (n_var < (k + 1)){
    stop("It takes at least k plus 1 variables to fit a k plus 1'th order t-cherry tree.")
  }

  tch_adj <- matrix(0, nrow = n_var, ncol = n_var)
  rownames(tch_adj) <- colnames(tch_adj) <- nodes
  for (i in 1:n_cliq) {
    tch_adj[tch_cliq[[i]], tch_cliq[[i]]] <- 1
    diag(tch_adj[tch_cliq[[i]], tch_cliq[[i]]]) <- 0
  }

  if (! all(gRbase::triangulateMAT(tch_adj) == tch_adj)){
    stop(paste("The cliques do not come from a triangulated graph.",
             "The cliques should correspond to a k'th order t-cherry",
             "tree so it must be triangulated.", collapse = " "))
  }

  if (sum(tch_adj) / 2 != (k - 1) * n_var - (1 / 2) * (k - 1) * k){
    stop(paste("The graph corresponding to the cliques does not have",
               "the correct number of edges for a k'th order t-cherry",
               "tree.", collapse = " "))
  }

  n_edges <- sum(tch_adj) / 2

  # Making all possible structures.

  first_cliques <- utils::combn(nodes, k + 1)
  models <- list()

  for (i in 1:ncol(first_cliques)){
    adj_matrix <- tch_adj
    adj_matrix[first_cliques[, i], first_cliques[, i]] <- 1
    diag(adj_matrix[first_cliques[, i], first_cliques[, i]]) <- 0

    if ((sum(adj_matrix) / 2 ) == (n_edges + 1)){
    model <- list("cliques" = list(first_cliques[, i]),
         "separators" = list(),
         "unused" = setdiff(nodes, first_cliques[, i]),
         "adj_matrix" = adj_matrix)
    models <- c(models, list(model))
    }
  }

  n_iter <- n_var - (k + 1)

  if (n_iter != 0){
    for (iter in 1:n_iter) {
      new_models <- models
      idx <- 1
      for (model in models) {
        for (clique in model$cliques) {
          seps <- utils::combn(clique, k)
          seps <- split(seps, rep(1:ncol(seps), each = nrow(seps)))
          for (sep in seps) {
            for (var in model$unused) {
              new_model <- list()
              new_clique <- c(sep, var)

              n_edges <- sum(model$adj_matrix) / 2

              new_matrix <- model$adj_matrix
              new_matrix[new_clique, new_clique] <- 1
              diag(new_matrix[new_clique, new_clique]) <- 0

              if (sum(new_matrix) / 2 == n_edges + 1){

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
