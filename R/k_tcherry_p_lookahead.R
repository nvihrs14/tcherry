#' Determine a k'th order t-cherry tree from data by adding p cliques
#' at a time.
#'
#' @description Determine the structure of a k'th order t-cherry tree
#' from data with realisations of n variables based on a greedy stepwise
#' approach. In each step p new cliques are added.
#'
#' @param data The data the tree structure should be based on.
#' @param k The order of the t-cherry tree.
#' @param p The number of new cliques considered in one step.
#' @param ... Additional arguments passed to \code{MIk}.
#'
#' @details Notice that for \eqn{p = 1} it is the same as using
#' \code{k_tcherry_step} and for \eqn{p = n - (k - 1)} it is the same
#' as using \code{tcherry_complete_search}.
#'
#' The algorithm for constructing the t-cherry tree from
#' data is based on an atempt to minimize the Kullback-Leibler
#' divergence by maximizing the weight
#' \deqn{\sum MI(clique) - \sum MI(separator).}
#' The first sum is over the cliques and the second over the
#' separators of the junction tree of the preliminary t-cherry tree.
#' In each step all possibilities of p new cliques are added to the
#' preliminary tree.
#' The one with the highest weight is chosen as the new preliminary
#' t-cherry tree, and the procedure is repeated untill all variables
#' has been added.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{adj_matrix} The adjacency matrix for the k'th order
#' t-cherry tree.
#' \item \code{weight} The weight of the final k'th order t-cherry tree.
#' \item \code{cliques} A list containing the cliques of
#'  the k'th order t-cherry tree.
#' \item \code{separators} A list containing the separators of a
#' junction tree for the k'th order t-cherry tree.
#' \item \code{n_edges} The number of edges in the resulting graph.
#' }
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @seealso \code{\link{MIk}} for mutual
#' information of k variables.
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
#'
#' data <- data.frame("var1" = as.character(var1),
#'                    "var2" = as.character(var2),
#'                    "var3" = as.character(var3),
#'                    "var4" = as.character(var4),
#'                    "var5" = as.character(var5),
#'                    "var6" = as.character(var6))
#'
#' # smooth used in MIk
#' (tch <- k_tcherry_p_lookahead(data, k = 3, p = 2, smooth = 0.1))
#'
#' # For plotting
#' library(gRbase)
#' library(Rgraphviz)
#' tcherry_tree <- as(tch$adj_matrix, "graphNEL")
#' plot(tcherry_tree)
#'
#' # For probability propagation
#' library(gRain)
#' model <- grain(tcherry_tree, data = data, smooth = 0.1)
#' querygrain(model)
#' @export

k_tcherry_p_lookahead <- function(data, k, p, ...){

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

  if (length(k) != 1 | length(p) != 1){
    stop("k and p must be single positive integers.")
  }

  if (p %% 1 != 0 | p < 1){
    stop("p must be a positive integer.")
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
    weight = MIk(col, data, ...)
    list("cliques" = list(col),
         "separators" = list(),
         "unused" = setdiff(nodes, col),
         "adj_matrix" = adj_matrix,
         "weight" = weight)
  })

  n_cliques_remaining <- n_var - k

  if (n_cliques_remaining == 0){
    weights <- sapply(models, function(l){
      l$weight
    })
    idx.max <- which.max(weights)
    model <- models[[idx.max]]
    model <- model[- 3]

    return(model)
  }

  first_time <- TRUE

  while (n_cliques_remaining > 0){

    if (first_time){
      n_iter <- min(p - 1, n_cliques_remaining)
    }else{
      n_iter <- min(p, n_cliques_remaining)
    }

    if (n_iter == 0){
      iterations <- c()
    }else{
      iterations <- 1:n_iter
    }

    for (iter in iterations) {
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

              new_model$cliques <- c(model$cliques,
                                         list(new_clique))
              new_model$new_cliques <- c(model$new_cliques,
                                         list(new_clique))
              new_model$separators <- c(model$separators,
                                            list(sep))
              new_model$new_separators <- c(model$new_separators,
                                        list(sep))
              new_model$unused <- setdiff(model$unused, var)
              new_model$adj_matrix <- new_matrix
              new_model$weight <- model$weight

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

    if (n_iter == 0){
      weights <- sapply(models, function(l){
      l$weight
      })
    }else{
      weights <- sapply(models, function(model){
      MIcliq <- sapply(model$new_cliques, MIk, data = data, ...)
      MIsep <- sapply(model$new_separators, MIk, data = data, ...)
      model$weight + sum(MIcliq) - sum(MIsep)
      })
    }

  idx.max <- which.max(weights)
  current_model <- models[[idx.max]]
  current_model$weight <- weights[idx.max]

  if (n_iter != 0){
    current_model <- current_model[-c(2, 4)]
  }
  models <- list(current_model)

  n_cliques_remaining <- n_cliques_remaining - n_iter
  first_time <- FALSE

  }

  current_model$n_edges <- sum(current_model$adj_matrix) / 2

  return(c(current_model[- 3]))

}
