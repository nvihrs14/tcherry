#' Determine the number of different edges for graphs.
#'
#' @description Determines the number of different edges for two undirected
#' graphs, with no loops, represented by adjacency matrices.
#'
#' @param adj_mat_1,adj_mat_2 Adjacency matrices for the two graphs.
#'
#' @return The number of edges where the two graphs differ.
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @examples
#' m1 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)
#' m2_same <- m1
#' m2_diff <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3, ncol = 3)
#'
#' diff_edges(m1, m2_same)
#' diff_edges(m1, m2_diff)
#'
#' @export

diff_edges <- function(adj_mat_1, adj_mat_2){

  if(! is.matrix(adj_mat_1) | ! is.matrix(adj_mat_2)){
    stop("Arguments must be matrices.")
  }

  if(! is.numeric(adj_mat_1) | ! is.numeric(adj_mat_2)){
    stop("Arguments must be numeric.")
  }

  if(any(!c(adj_mat_1, adj_mat_2) %in% 0:1)){
    stop("Arguments must be adjacency matrices for unweighted graphs.
         Therefore all entries must be 0 or 1.")
  }

  if(! isSymmetric(adj_mat_1) | ! isSymmetric(adj_mat_2)){
    stop("Only undirected graphs are supported so arguments must be
         symmetric.")
  }

  if(any(c(diag(adj_mat_1), diag(adj_mat_2)) != 0)){
    stop("Loops are not supported so diagonal must be all 0.")
  }

  u_tri_1_idx <- upper.tri(adj_mat_1)
  u_tri_2_idx <- upper.tri(adj_mat_2)

  u_tri_1 <- adj_mat_1[u_tri_1_idx]
  u_tri_2 <- adj_mat_2[u_tri_2_idx]

  is_different <- u_tri_1 != u_tri_2
  length(which(is_different))
}
