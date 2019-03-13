tcherry_complete_search <- function(data, k, ...){
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
      new_models <- models
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

  return(list("model" = model,
              "n_models" = length(models)))

}

