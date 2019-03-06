#' Determine a t-cherry tree from data
#'
#' @description Determine the structure of a t-cherry tree
#' from data based on a greedy stepwise approach.
#'
#' @param data The data the tree structure should be based on.
#' @param ... Additional arguments passed to \code{MI2} and
#' \code{MI3}.
#'
#' @details The algorithm for constructing the t-cherry tree from
#' data is based on an atempt to minimize the Kullback-Leibler
#' divergence. The first cherry is chosen as the triplet with
#' highest mutual information. This is the preliminary t-cherry
#' tree. Then all possible new cherries are added stepwise to this
#' tree and the weight \deqn{\sum MI3(Clique) - \sum MI2(Separator)}
#' is calculated.
#' The first sum is over the cliques and the second over the
#' separators of the junction tree of the preliminary t-cherry tree.
#' The one with the highest weight is chosen as the new preliminary
#' t-cherry tree, and the procedure is repeated untill all variables
#' has been added.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{adj_matrix} The adjacency matrix for the t-cherry
#' tree.
#' \item \code{weight} The weight of the final t-cherry tree.
#' \item \code{cliques} A list containing the cliques (cherries) of
#'  the t-cherry tree.
#' \item \code{separators} A list containing the separators of a
#' junction tree for the t-cherry tree.
#' }
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @seealso \code{\link{MI2}} and \code{\link{MI3}} for mutual
#' information of two and three variables respectively.
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
#' # smooth used in both MI2 and MI3
#' (tch <- tcherry_step(data, smooth = 0.1))
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

tcherry_step <- function(data, ...){
  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix.")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors.")
  }

  nodes <- colnames(data)
  n_var <- length(nodes)
  adj_matrix <- matrix(0, nrow = n_var, ncol = n_var,
                       dimnames = list(nodes, nodes))
  cliques <- as.list(rep(NA, n_var - 2))
  separators <- as.list(rep(NA, n_var - 3))
  # For use in for-loop
  MI2_fun <- function(var1, var2){
    MI2(data[, var1], data[, var2], ...)
  }

  # Adding first cherry
  triples <- t(utils::combn(nodes, 3))

  MI3_fun <- function(var1, var2, var3){
    MI3(data[, var1], data[, var2], data[, var3], ...)
  }

  MI3 <- mapply(MI3_fun, triples[, 1], triples[, 2], triples[, 3])
  MI3_tab <- data.frame(triples, MI3, stringsAsFactors = FALSE)

  idx.max <- which.max(MI3_tab$MI3)

  tcherry_nodes <- unlist(MI3_tab[idx.max, 1:3])
  names(tcherry_nodes) <- NULL
  cliques[[1]] <- tcherry_nodes
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

  weight <- MI3_tab$MI3[idx.max]

  idx_list <- 1

  while (length(tcherry_nodes) != n_var) {
    weight_next_step <- rep(NA, nrow(tcherry_edges) *
                             length(nodes_remaining))
    new_cliques_list <- new_seps_list <- new_var_list <-
      as.list(weight_next_step)
    idx <- 1
    for (i in 1:nrow(tcherry_edges)){
      for (var in nodes_remaining){
        new_sep <- unlist(tcherry_edges[i, ])
        names(new_sep) <- NULL
        new_clique <- unlist(c(new_sep, var))
        MI2_sep <- MI2_fun(new_sep[1], new_sep[2])
        MI3_clique <- MI3_fun(new_clique[1], new_clique[2],
                              new_clique[3])
        weight_next_step[idx] <- weight + MI3_clique - MI2_sep
        new_cliques_list[[idx]] <- new_clique
        new_seps_list[[idx]] <- new_sep
        new_var_list[[idx]] <- var
        idx <- idx + 1
      }
    }

    idx_max_weight <- which.max(weight_next_step)

    weight <- weight_next_step[idx_max_weight]

    tcherry_nodes <- unlist(c(tcherry_nodes,
                              new_var_list[[idx_max_weight]]))
    nodes_remaining <- nodes[! nodes %in% tcherry_nodes]

    cliques[[idx_list + 1]] <- new_cliques_list[[idx_max_weight]]
    separators[[idx_list]] <- new_seps_list[[idx_max_weight]]
    idx_list <- idx_list + 1

    tcherry_edges[nrow(tcherry_edges) + 1, ] <-
      c(new_var_list[[idx_max_weight]],
        new_seps_list[[idx_max_weight]][1])
    tcherry_edges[nrow(tcherry_edges) + 1, ] <-
      c(new_var_list[[idx_max_weight]],
        new_seps_list[[idx_max_weight]][2])

    edge_1 <- tcherry_edges[nrow(tcherry_edges) - 1, 1]
    edge_2 <- tcherry_edges[nrow(tcherry_edges) - 1, 2]
    edge_3 <- tcherry_edges[nrow(tcherry_edges), 2]
    adj_matrix[edge_1, edge_2] <-
      adj_matrix[edge_2, edge_1] <-
      adj_matrix[edge_1, edge_3] <-
      adj_matrix[edge_3, edge_1] <- 1
  }
  return(list("adj_matrix" = adj_matrix,
              "weight" = weight,
              "cliques" = cliques,
              "separators" = separators))
}
