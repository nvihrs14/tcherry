k_tcherry_step <- function(data, k, ...){
  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix.")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors.")
  }

  if (k %% 1 != 0 | k <= 1){
    stop("k must be a positive integer and at least 2.")
  }

  nodes <- colnames(data)
  n_var <- length(nodes)
  adj_matrix <- matrix(0, nrow = n_var, ncol = n_var,
                       dimnames = list(nodes, nodes))
  cliques <- as.list(rep(NA, n_var - (k - 1)))
  separators <- as.list(rep(NA, n_var - k))

  # Adding first cherry
  poss_cliq <- utils::combn(nodes, k)
  poss_cliq <- split(poss_cliq, rep(1:ncol(poss_cliq),
                                    each = nrow(poss_cliq)))

  MI <- sapply(poss_cliq, MIk, data = data, ...)
  idx.max <- which.max(MI)

  tcherry_nodes <- poss_cliq[[idx.max]]
  cliques[[1]] <- tcherry_nodes
  nodes_remaining <- setdiff(nodes, tcherry_nodes)
  tcherry_hyperedges <- utils::combn(tcherry_nodes, k-1)
  tcherry_hyperedges <- split(tcherry_hyperedges,
                              rep(1:ncol(tcherry_hyperedges),
                                  each = nrow(tcherry_hyperedges)))

  adj_matrix[tcherry_nodes, tcherry_nodes] <- 1
  diag(adj_matrix[tcherry_nodes, tcherry_nodes]) <- 0

  score <- max(MI)

  idx_list <- 1

  while (length(tcherry_nodes) != n_var) {
    score_next_step <- rep(NA, length(tcherry_hyperedges) *
                             length(nodes_remaining))
    new_cliques_list <- new_seps_list <- new_var_list <-
      as.list(score_next_step)
    idx <- 1
    for (i in 1:length(tcherry_hyperedges)){
      for (var in nodes_remaining){
        new_sep <- tcherry_hyperedges[[i]]
        new_clique <- c(new_sep, var)

        MI_sep <- MIk(new_sep, data, ...)
        MI_clique <- MIk(new_clique, data, ...)

        score_next_step[idx] <- score + MI_clique - MI_sep
        new_cliques_list[[idx]] <- new_clique
        new_seps_list[[idx]] <- new_sep
        new_var_list[[idx]] <- var
        idx <- idx + 1
      }
    }

    idx_max_score <- which.max(score_next_step)
    score <- score_next_step[idx_max_score]

    new_clique <- new_cliques_list[[idx_max_score]]
    new_sep <- new_seps_list[[idx_max_score]]
    new_var <- new_var_list[[idx_max_score]]

    tcherry_nodes <- c(tcherry_nodes, new_var)
    nodes_remaining <- setdiff(nodes, tcherry_nodes)

    cliques[[idx_list + 1]] <- new_clique
    separators[[idx_list]] <- new_sep
    idx_list <- idx_list + 1

    new_hyper_edges <- utils::combn(new_clique, k - 1)
    new_hyper_edges <- split(new_hyper_edges,
                                rep(1:ncol(new_hyper_edges),
                                    each = nrow(new_hyper_edges)))
    tcherry_hyperedges <- c(tcherry_hyperedges, new_hyper_edges)
    tcherry_hyperedges <- unique(tcherry_hyperedges)

    adj_matrix[new_clique, new_clique] <- 1
    diag(adj_matrix[new_clique, new_clique]) <- 0

  }

  return(list("adj_matrix" = adj_matrix,
              "score" = score,
              "cliques" = cliques,
              "separators" = separators))
}
