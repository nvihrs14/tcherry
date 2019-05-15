#' @rdname increase_order1
#' @export

increase_order2 <- function(tch_cliq, data, ...){

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

  nodes <- colnames(data)
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

  cliques <- as.list(rep(NA, n_var - k))
  separators <- as.list(rep(NA, n_var - (k + 1)))

  # Adding first cherry.
  poss_cliq <- utils::combn(nodes, k + 1)
  poss_cliq <- split(poss_cliq, rep(1:ncol(poss_cliq),
                                    each = nrow(poss_cliq)))
  dat_first_cherry <- data.frame(cliq = I(poss_cliq),
                                 MI = rep(NA, length(poss_cliq)),
                                 adj_mat = I(as.list(
                                   rep(NA, length(poss_cliq)))))
  n_edges <- sum(tch_adj) / 2

  for (i in 1:nrow(dat_first_cherry)) {
    adj_matrix_temp <- tch_adj
    adj_matrix_temp[dat_first_cherry$cliq[[i]],
                    dat_first_cherry$cliq[[i]]] <- 1
    diag(adj_matrix_temp[dat_first_cherry$cliq[[i]],
                         dat_first_cherry$cliq[[i]]]) <- 0
    if ((sum(adj_matrix_temp) / 2 ) == (n_edges + 1)){
      dat_first_cherry$MI[i] <- MIk(dat_first_cherry$cliq[[i]],
                                      data, ...)
      dat_first_cherry$adj_mat[[i]] <- adj_matrix_temp
    }
  }

  idx.max <- which.max(dat_first_cherry$MI)

  tcherry_nodes <- dat_first_cherry$cliq[[idx.max]]
  cliques[[1]] <- tcherry_nodes
  nodes_remaining <- setdiff(nodes, tcherry_nodes)

  tcherry_hyperedges <- utils::combn(tcherry_nodes, k)
  tcherry_hyperedges <- split(tcherry_hyperedges,
                              rep(1:ncol(tcherry_hyperedges),
                                  each = nrow(tcherry_hyperedges)))

  adj_matrix <- dat_first_cherry$adj_mat[[idx.max]]

  weight <- dat_first_cherry$MI[idx.max]

  n_edges <- n_edges + 1

  # Adding remaining cherries.

  idx.dat <- 1
  idx.list <- 1

  n_nodes_remaining_median <- floor((n_var - (k + 1)) / 2 + 1)
  n_hyp_edges_median <- 1 + (k) * (
    (n_var - n_nodes_remaining_median) - (k))

  weight_cliq_sep <- MI_cliq <- MI_sep <- new_var <-
    rep(NA, n_nodes_remaining_median * n_hyp_edges_median)

  new_cliques_list <- new_seps_list <-
    as.list(weight_cliq_sep)

  dat_new_poss <- data.frame(new_cliq = I(new_cliques_list),
                             new_sep = I(new_seps_list),
                             new_var = new_var,
                             MI_cliq = MI_cliq,
                             MI_sep = MI_sep,
                             weight_increase = weight_cliq_sep)

  while (length(tcherry_nodes) != n_var){
    for (i in 1:length(tcherry_hyperedges)){
      for (var in nodes_remaining){
        new_sep <- tcherry_hyperedges[[i]]
        dat_new_poss$new_sep[[idx.dat]] <- new_sep
        new_cliq <- c(dat_new_poss$new_sep[[idx.dat]], var)
        dat_new_poss$new_cliq[[idx.dat]] <- new_cliq

        adj_matrix_temp <- adj_matrix
        adj_matrix_temp[new_cliq, new_cliq] <- 1
        diag(adj_matrix_temp[new_cliq, new_cliq]) <- 0

        if ((sum(adj_matrix_temp) / 2 ) == (n_edges + 1)){
          dat_new_poss$MI_sep[idx.dat] <- MIk(new_sep, data, ...)
          dat_new_poss$MI_cliq[idx.dat] <- MIk(new_cliq, data, ...)

          dat_new_poss$weight_increase[idx.dat] <-
          dat_new_poss$MI_cliq[idx.dat] - dat_new_poss$MI_sep[idx.dat]
        }

        dat_new_poss$new_var[[idx.dat]] <- var
        idx.dat <- idx.dat + 1
      }
    }

    idx_max_weight <- which.max(dat_new_poss$weight_increase)
    weight <- weight + dat_new_poss$weight_increase[idx_max_weight]

    new_clique <- dat_new_poss$new_cliq[[idx_max_weight]]
    new_sep <- dat_new_poss$new_sep[[idx_max_weight]]
    new_var <- dat_new_poss$new_var[idx_max_weight]

    tcherry_nodes <- c(tcherry_nodes, new_var)
    nodes_remaining <- setdiff(nodes, tcherry_nodes)

    cliques[[idx.list + 1]] <- new_clique
    separators[[idx.list]] <- new_sep
    idx.list <- idx.list + 1

    new_hyper_edges <- utils::combn(new_clique, k)
    new_hyper_edges <- split(new_hyper_edges,
                             rep(1:ncol(new_hyper_edges),
                                 each = nrow(new_hyper_edges)))
    idx.new <- sapply(new_hyper_edges, function(e){
      length(setdiff(e, new_var)) != k
    })
    tcherry_hyperedges <- new_hyper_edges[idx.new]

    adj_matrix[new_clique, new_clique] <- 1
    diag(adj_matrix[new_clique, new_clique]) <- 0

    idx.new.var <- dat_new_poss$new_var == new_var &
      !is.na(dat_new_poss$new_var)
    dat_new_poss[idx.new.var, ] <- NA
    ord <- order(dat_new_poss$new_var)
    dat_new_poss <- dat_new_poss[ord, ]

    idx.dat <- which(is.na(dat_new_poss$new_var))[1]
    n_edges <- n_edges + 1
  }

  n_edges_graph <- sum(adj_matrix) / 2

  return(list("adj_matrix" = adj_matrix,
              "weight" = weight,
              "cliques" = cliques,
              "separators" = separators,
              "n_edges" = n_edges_graph))
}
