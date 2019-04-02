edges_in_one_clique <- function(cliques){
  edges <- lapply(cliques, function(cliq){
    if (length(cliq) > 1){
      comb <- utils::combn(cliq, 2)
      split(comb, rep(1:ncol(comb), each = nrow(comb)))
    }
  })

  n_edges_cliq <- sapply(edges, length)
  edges <- unlist(edges, recursive = FALSE)
  n_edges <- length(unique(edges))

  data_edges <- data.frame(edge = I(edges),
                           clique = I(rep(cliques, n_edges_cliq)))


  not_in_one_cliq <- unique(edges[duplicated(edges)])
  data_edges_yes <- data_edges[! data_edges$edge %in% not_in_one_cliq, ]
  list("data" = data_edges_yes,
       "n" = n_edges)
}
