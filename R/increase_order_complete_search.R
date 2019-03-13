increase_order_complete_search <- function(tch_cliq, data, ...){
  nodes <- names(data)
  n_var <- length(nodes)

  n_cliq <- length(tch_cliq)
  k <- length(tch_cliq[[1]])

  tch_adj <- matrix(0, nrow = n_var, ncol = n_var)
  rownames(tch_adj) <- colnames(tch_adj) <- nodes
  for (i in 1:n_cliq) {
    tch_adj[tch_cliq[[i]], tch_cliq[[i]]] <- 1
    diag(tch_adj[tch_cliq[[i]], tch_cliq[[i]]]) <- 0
  }

  n_edges <- sum(tch_adj) / 2

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

  return(list("model" = model,
              "n_models" = length(models)))

}
