
tcherry <- function(data, ...){
  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors")
  }

  CL <- ChowLiu(data, ...)
  tree <- CL$skeleton_adj

  nodes <- names(data)
  n_var <- length(nodes)
  triples <- t(combn(nodes, 3))

  MI3_fun <- function(var1, var2, var3){
    MI3(data[, var1], data[, var2], data[, var3], ...)
  }

  MI3 <- mapply(MI3_fun, triples[, 1], triples[, 2], triples[, 3])
  MI3_tab <- data.frame(triples, MI3, stringsAsFactors = FALSE)

  ord_idx <- order(MI3_tab$MI3, decreasing = TRUE)
  MI3_tab <- MI3_tab[ord_idx, ]
  rownames(MI3_tab) <- NULL

  MI <- MI3_tab

  n_edges <- sum(tree) / 2
  tcherry_nodes <- c()

  i <- 1
  while (length(tcherry_nodes) == 0) {
    adj_matrix_temp <- tree
    edge_1 <- MI3_tab[i,1]
    edge_2 <- MI3_tab[i,2]
    edge_3 <- MI3_tab[i,3]
    adj_matrix_temp[edge_1, edge_2] <-
      adj_matrix_temp[edge_2, edge_1] <-
      adj_matrix_temp[edge_1, edge_3] <-
      adj_matrix_temp[edge_3, edge_1] <-
      adj_matrix_temp[edge_2, edge_3] <-
      adj_matrix_temp[edge_3, edge_2] <- 1
    if ((sum(adj_matrix_temp) / 2 ) == (n_edges + 1)){
      adj_matrix <- adj_matrix_temp
      tcherry_nodes <- unlist(c(tcherry_nodes, MI3_tab[i, 1:3]))
      n_edges <- n_edges + 1
      MI3_tab <- MI3_tab[- i, ]
      cliques <- list(tcherry_nodes)
      }
    }

  i <- k <- 1
  j <- 2

  separators <- list()

  while (length(tcherry_nodes) < n_var) {
    n_var_in_tcherry <- length(which(MI3_tab[i, 1:3]
                                     %in% tcherry_nodes))
    if (n_var_in_tcherry == 2){
      adj_matrix_temp <- adj_matrix
      edge_1 <- MI3_tab[i,1]
      edge_2 <- MI3_tab[i,2]
      edge_3 <- MI3_tab[i,3]
      adj_matrix_temp[edge_1, edge_2] <-
        adj_matrix_temp[edge_2, edge_1] <-
        adj_matrix_temp[edge_1, edge_3] <-
        adj_matrix_temp[edge_3, edge_1] <-
        adj_matrix_temp[edge_2, edge_3] <-
        adj_matrix_temp[edge_3, edge_2] <- 1
      if ((sum(adj_matrix_temp) / 2 ) == (n_edges + 1)){
        adj_matrix <- adj_matrix_temp
        s_idx <- which(MI3_tab[i, 1:3] %in% tcherry_nodes)
        cliques[[j]] <- unlist(c(MI3_tab[i, 1:3]))
        separators[[k]] <- cliques[[j]][s_idx]


        tcherry_nodes <- unique(unlist(c(tcherry_nodes,
                                        MI3_tab[i, 1:3])))
        n_edges <- n_edges + 1


        idx_delete <- which(rowSums(matrix(
          as.matrix(MI3_tab[, 1:3]) %in% tcherry_nodes,
          nrow = nrow(MI3_tab))) == 3)

        MI3_tab <- MI3_tab[-idx_delete, ]

        i <- 0
        j <- j + 1
        k <- k + 1
        }
      }
    i <- i + 1
  }

  return(list("adj_tcherry" = adj_matrix,
              "cliques" = cliques,
              "separators" = separators))

}
