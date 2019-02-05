#' Determine t-cherry tree from a Chow-Liu tree
#'
#' @description Determine the structure of a t-cherry tree
#' from a Chow-Liu tree for data.
#'
#' @param data The data the tree structure should be based on.
#' @param ... Additional arguments passed to \code{MI2} and
#' \code{MI3}.
#'
#' @details The algorithm for constructing the t-cherry tree from
#' a Chow-Liu tree is as described in \insertRef{EKTS}{tcherry}.
#'
#' The algorithm is greedy in the sence that it always attempts to
#' use the three variables with highest mutual information as the
#' next cherry.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{adj_tcherry} The adjacency matrix for the t-cherry
#' tree.
#' \item \code{cliques} A list containing the cliques (cherries) of
#'  the t-cherry tree.
#' \item \code{separators} A list containing the separators of a
#' junction tree for the t-cherry tree.
#' \item \code{MI} A data frame containing the mutual information
#' of all triplets.
#' }
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @references
#' \insertRef{EKTS}{tcherry}
#'
#' @seealso \code{\link{ChowLiu}} for construction of Chow-Liu
#' trees and \code{\link{MI3}} for mutual information of three
#' variables.
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
#' (tch <- tcherry_CL(data, smooth = 0.1))
#'
#' # For plotting
#' library(gRbase)
#' library(Rgraphviz)
#' tcherry_tree <- as(tch$adj_tcherry, "graphNEL")
#' plot(tcherry_tree)
#'
#' # For probability propagation
#' library(gRain)
#' model <- grain(tcherry_tree, data = data, smooth = 0.1)
#' querygrain(model)
#' @export

tcherry_CL <- function(data, ...){
  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors")
  }

  CL <- ChowLiu(data, CPTs = FALSE, ...)
  tree <- CL$skeleton_adj

  nodes <- colnames(data)
  n_var <- length(nodes)
  triples <- t(utils::combn(nodes, 3))

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

  # Find the first cherry
  i <- 1
  while (length(tcherry_nodes) == 0) {
    adj_matrix_temp <- tree
    edge_1 <- MI3_tab[i, 1]
    edge_2 <- MI3_tab[i, 2]
    edge_3 <- MI3_tab[i, 3]
    adj_matrix_temp[edge_1, edge_2] <-
      adj_matrix_temp[edge_2, edge_1] <-
      adj_matrix_temp[edge_1, edge_3] <-
      adj_matrix_temp[edge_3, edge_1] <-
      adj_matrix_temp[edge_2, edge_3] <-
      adj_matrix_temp[edge_3, edge_2] <- 1
    if ( (sum(adj_matrix_temp) / 2 ) == (n_edges + 1)){
      adj_matrix <- adj_matrix_temp
      tcherry_nodes <- unlist(c(tcherry_nodes, MI3_tab[i, 1:3]))
      n_edges <- n_edges + 1
      MI3_tab <- MI3_tab[- i, ]
      cliques <- list(tcherry_nodes)
      names(cliques[[1]]) <- NULL
      }
    i <- i + 1
    }

  # Add remaining nodes via new cherries.
  i <- k <- 1
  j <- 2

  separators <- list()

  while (length(tcherry_nodes) < n_var) {
    n_var_in_tcherry <- length(which(MI3_tab[i, 1:3]
                                     %in% tcherry_nodes))
    if (n_var_in_tcherry == 2){
      adj_matrix_temp <- adj_matrix
      edge_1 <- MI3_tab[i, 1]
      edge_2 <- MI3_tab[i, 2]
      edge_3 <- MI3_tab[i, 3]
      adj_matrix_temp[edge_1, edge_2] <-
        adj_matrix_temp[edge_2, edge_1] <-
        adj_matrix_temp[edge_1, edge_3] <-
        adj_matrix_temp[edge_3, edge_1] <-
        adj_matrix_temp[edge_2, edge_3] <-
        adj_matrix_temp[edge_3, edge_2] <- 1
      if ( (sum(adj_matrix_temp) / 2 ) == (n_edges + 1)){
        adj_matrix <- adj_matrix_temp
        s_idx <- which(MI3_tab[i, 1:3] %in% tcherry_nodes)
        cliques[[j]] <- unlist(c(MI3_tab[i, 1:3]))
        separators[[k]] <- cliques[[j]][s_idx]
        names(cliques[[j]]) <- names(separators[[k]]) <- NULL

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
              "separators" = separators,
              "MI" = MI))

}
