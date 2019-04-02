
thinning_edges <- function(cliques, data, alpha = 0.05, ...){

  edge_information <- edges_in_one_clique(cliques)
  data_edges_yes <- edge_information$data

  n_edges <- edge_information$n

  have_tested <- as.list(rep(NA), n_edges)
  n_edges_removed <- 0

  i <- 1
  j <- 1

  while (i <= nrow(data_edges_yes)){
    edge <- data_edges_yes$edge[[i]]
    clique <- data_edges_yes$clique[[i]]

    var1 <- edge[1]
    var2 <- edge[2]
    cond <- setdiff(clique, edge)

    test <- cond_independence_test(var1, var2, cond, data, ...)

    if (test$p_value < alpha){
      have_tested[[j]] <- edge
      i <- i + 1
      j <- j + 1
    }

    if (test$p_value >= alpha){
      n_edges_removed <- n_edges_removed + 1

      cliques <- setdiff(cliques, list(clique))
      new_cliques <- list()
      new_clique_1 <- c(var1, cond)

      is_subset_1 <- lapply(cliques, function(cliq){
        length(setdiff(new_clique_1, cliq)) == 0
      })

      if (! any(unlist(is_subset_1))){
        new_cliques <- c(new_cliques, list(new_clique_1))
      }

      new_clique_2 <- c(var2, cond)

      is_subset_2 <- lapply(cliques, function(cliq){
        length(setdiff(new_clique_2, cliq)) == 0
      })

      if (! any(unlist(is_subset_2))){
        new_cliques <- c(new_cliques, list(new_clique_2))
      }

      cliques <- c(cliques, new_cliques)

      data_edges_yes <- edges_in_one_clique(cliques)$data
      data_edges_yes <-
        data_edges_yes[! data_edges_yes$edge %in% have_tested, ]

      i <- 1
    }
  }

  # Reconstructing adjacency matrix from cliques

  nodes <- sort(unique(unlist(cliques)))
  n_var <- length(nodes)
  n_cliq <- length(cliques)

  adj_matrix <- matrix(0, nrow = n_var, ncol = n_var)
  rownames(adj_matrix) <- colnames(adj_matrix) <- nodes
  for (i in 1:n_cliq) {
    if (length(cliques[[i]]) != 1){
    adj_matrix[cliques[[i]], cliques[[i]]] <- 1
    diag(adj_matrix[cliques[[i]], cliques[[i]]]) <- 0
    }
  }

  return(list("adj_matrix" = adj_matrix,
              "cliques" = cliques,
              "n_edges_removed" = n_edges_removed))
}






cliques <- list(c("var1", "var2", "var3"),
                c("var2", "var3", "var5"),
                c("var5", "var6", "var7"),
                c("var1", "var4"),
                c("var2", "var5", "var6"))
