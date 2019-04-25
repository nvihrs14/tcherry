#' Thinning of edges in a graphical model with a triangulated graph
#'
#' @description Thinning of edges in a graphical model with a
#' triangulated graph based on a likelihood ratio test for conditional
#' independence.
#'
#' @param cliques A list containing the cliques of the triangulated
#' graph, which should be thinned.
#' @param separators A list containing the separators of the junction
#' tree of the triangulated graph, which should be thinned.
#' @param data Data with realisations of the variables in the graph.
#' @param alpha Significance level used in the conditional
#' independece test.
#' @param ... Additional arguments passed to
#' \code{cond_independence_test}
#'
#' @details
#' The edges in the graph are removed one by one if the respetive
#' conditional independence test cannot be rejected. An edge is only
#' considered if its removal results in a new triangulated graph. This
#' means that only edges in one clique only are considered.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{adj_matrix} The adjacency matrix of the thinned graph.
#' \item \code{cliques} The cliques of the thinned graph.
#' \item \code{separators} The separators of the junction tree for
#' the thinned graph.
#' \item \code{n_edges} The number of edges in the resulting graph.
#' \item \code{n_edges_removed} The number of removed edges.
#' }
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @seealso \code{\link{cond_independence_test}} for the test used.
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
#' cliques <- list(c("var1", "var2", "var3"),
#'                 c("var2", "var3", "var5"),
#'                 c("var5", "var6", "var7"),
#'                 c("var1", "var4"),
#'                 c("var2", "var5", "var6"))
#'
#' separators <- list(c("var1"),
#'                    c("var2", "var3"),
#'                    c("var2", "var5"),
#'                    c("var5", "var6"))
#'
#' thinning_edges(cliques, separators, data = data, alpha = 0.1,
#'                smooth = 0.1)
#'
#' @export

thinning_edges <- function(cliques, separators, data, alpha = 0.05, ...){

  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix.")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors.")
  }

  if (! is.list(cliques)){
    stop(paste("Cliques must be given in a list, each entry containing",
               "a vector with the names of the variables in the clique.",
               collapse = " "))
  }

  if (! compare::compare(unique(unlist(cliques)), colnames(data),
                         ignoreOrder = TRUE)$result){
    stop(paste("The column names of data must be the same as the",
               "variable names in tch_cliq. All variables in data must",
               "be in at least one clique.", collapse = " "))
  }

  if (! is.list(separators)){
  stop(paste("Separators must be given in a list, each entry containing",
           "a vector with the names of the variables in the separator.",
               collapse = " "))
  }

  if (length(setdiff(unique(unlist(separators)), colnames(data))) != 0){
    stop("All variable names in separators should be in data.")
  }

  if (length(alpha) > 1){
    stop("alpha must be a single non-negative value.")
  }
  else if (!is.numeric(alpha)) {
    stop("alpha must be numeric.")
  }
  else if (alpha <= 0){
    stop("alpha must be a positive numeric value.")
  }

  cliques <- lapply(cliques, sort)
  separators <- lapply(separators, sort)

  edge_information <- edges_in_one_clique(cliques)
  data_edges_yes <- edge_information$data
  data_edges_yes$edge <- lapply(data_edges_yes$edge, sort)

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
      new_clique_1 <- sort(c(var1, cond))

      is_sep_1 <- sapply(separators, function(sep){
        compare::compare(new_clique_1, sep)$result
      })

      if (! any(is_sep_1)){
        new_cliques <- c(new_cliques, list(new_clique_1))
      } else {
        idx_sep <- which(is_sep_1)[1]
        separators <- separators[- idx_sep]
      }


      new_clique_2 <- sort(c(var2, cond))

      is_sep_2 <- sapply(separators, function(sep){
        compare::compare(new_clique_2, sep)$result
      })

      if (! any(is_sep_2)){
        new_cliques <- c(new_cliques, list(new_clique_2))
      } else {
        idx_sep <- which(is_sep_2)[1]
        separators <- separators[- idx_sep]
      }

      cliques <- c(cliques, new_cliques)
      new_sep <- intersect(new_clique_1, new_clique_2)
      separators <- c(separators, list(new_sep))

      data_edges_yes <- edges_in_one_clique(cliques)$data
      data_edges_yes$edge <- lapply(data_edges_yes$edge, sort)
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

  n_edges_graph <- sum(adj_matrix) / 2

  return(list("adj_matrix" = adj_matrix,
              "cliques" = cliques,
              "separators" = separators,
              "n_edges" = n_edges_graph,
              "n_edges_removed" = n_edges_removed))
}







