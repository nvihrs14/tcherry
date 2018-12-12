
tcherry_step <- function(data, ...){
  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors")
  }

  nodes <- colnames(data)
  n_var <- length(nodes)
  adj_matrix <- matrix(0, nrow = n_var, ncol = n_var,
                       dimnames = list(nodes, nodes))

  # MI for all pairs
  pairs <- t(utils::combn(nodes, 2))

  MI2_fun <- function(var1, var2){
    MI2(data[, var1], data[, var2], ...)
  }

  MI2 <- mapply(MI2_fun, pairs[, 1], pairs[, 2])
  MI2_tab <- data.frame(pairs, MI2, stringsAsFactors = FALSE)

  # Adding first cherry
  triples <- t(utils::combn(nodes, 3))

  MI3_fun <- function(var1, var2, var3){
    MI3(data[, var1], data[, var2], data[, var3], ...)
  }

  MI3 <- mapply(MI3_fun, triples[, 1], triples[, 2], triples[, 3])
  MI3_tab <- data.frame(triples, MI3, stringsAsFactors = FALSE)

  idx.max <- which.max(MI3_tab$MI3)

  tcherry_nodes <- MI3_tab[idx.max, 1:3]
  nodes_remaining <- nodes[! nodes %in% tcherry_nodes]
  edge_1 <- MI3_tab[idx.max, 1]
  edge_2 <- MI3_tab[idx.max, 2]
  edge_3 <- MI3_tab[idx.max, 3]
  e1 <- c(edge_1, edge_1, edge_2)
  e2 <- c(edge_2, edge_3, edge_3)
  tcherry_edges <- data.frame(e1, e2, stringsAsFactors = FALSE)
  adj_matrix[edge_1, edge_2] <-
    adj_matrix[edge_2, edge_1] <-
    adj_matrix[edge_1, edge_3] <-
    adj_matrix[edge_3, edge_1] <-
    adj_matrix[edge_2, edge_3] <-
    adj_matrix[edge_3, edge_2] <- 1

  score <- MI3_tab$MI3[idx.max]

  while (length(tcherry_nodes) != n_var) {
    score_next_step <- rep(NA, nrow(tcherry_edges) *
                             length(nodes_remaining))
    new_cliques_list <- new_seps_list <- new_var_list <-
      as.list(score_next_step)
    idx <- 1
    for (i in 1:nrow(tcherry_edges)){
      for (var in nodes_remaining){
        new_sep <- tcherry_edges[i, ]
        new_clique <- unlist(c(new_sep, var))
        sep_in_row <- apply(MI2_tab, 1, function(r){
          new_sep[1] %in% r & new_sep[2] %in% r
        })
        MI2_sep <- MI2_tab$MI2[sep_in_row]
        clique_in_row <- apply(MI3_tab, 1, function(r){
          new_clique[1] %in% r & new_clique[2] %in% r &
            new_clique[3] %in% r
        })
        MI3_clique <- MI3_tab$MI3[clique_in_row]
        score_next_step[idx] <- score + MI3_clique - MI2_sep
        new_cliques_list[[idx]] <- new_clique
        new_seps_list[[idx]] <- new_sep
        new_var_list[[idx]] <- var
        idx <- idx + 1
      }
    }
    idx_max_score <- which.max(score_next_step)
    score <- score_next_step[idx_max_score]
    tcherry_nodes <- c(tcherry_nodes, new_var_list[[idx_max_score]])
    nodes_remaining <- nodes[! nodes %in% tcherry_nodes]
    tcherry_edges[nrow(tcherry_edges) + 1, ] <-
      c(new_var_list[[idx_max_score]],
        new_seps_list[[idx_max_score]][1])
    tcherry_edges[nrow(tcherry_edges) + 1, ] <-
      c(new_var_list[[idx_max_score]],
        new_seps_list[[idx_max_score]][2])
    edge_1 <- tcherry_edges[nrow(tcherry_edges) - 1, 1]
    edge_2 <- tcherry_edges[nrow(tcherry_edges) - 1, 2]
    edge_3 <- tcherry_edges[nrow(tcherry_edges), 2]
    adj_matrix[edge_1, edge_2] <-
      adj_matrix[edge_2, edge_1] <-
      adj_matrix[edge_1, edge_3] <-
      adj_matrix[edge_3, edge_1] <- 1
  }
  return(list("adj_matrix" = adj_matrix,
              "score" = score))
}


