#' Determine the number of differing edges for t-cherry trees
#'
#' @description Determines the number of differing edges for two t-cherry
#' trees over the same universe and of the same order represented by
#' adjacency matrices,
#' i.e. the number of edges in each graph minus the number of edges
#' common to both t-cherry trees.
#'
#' @param adj_mat_1,adj_mat_2 Adjacency matrices for the two graphs.
#'
#' @return The number of edges for which the two graphs differ.
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @examples
#' m1 <- matrix(c(0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0),
#'              nrow = 4, ncol = 4)
#' m2 <- matrix(c(0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0),
#'              nrow = 4, ncol = 4)
#' colnames(m1) <- rownames(m1) <- letters[1:4]
#' colnames(m2) <- rownames(m2) <- letters[1:4]
#'
#' diff_edges_tch(m1, m1)
#' diff_edges_tch(m1, m2)
#' @export

diff_edges_tch <- function(adj_mat_1, adj_mat_2){

  if (! is.matrix(adj_mat_1) | ! is.matrix(adj_mat_2)){
    stop("Arguments must be matrices.")
  }

  if (! isTRUE(all.equal(dim(adj_mat_1), dim(adj_mat_2)))){
    stop("The matrices must have the same dimensions.")
  }

  if (! is.numeric(adj_mat_1) | ! is.numeric(adj_mat_2)){
    stop("Arguments must be numeric.")
  }

  if (any(! c(adj_mat_1, adj_mat_2) %in% 0:1)){
    stop("Arguments must be adjacency matrices for unweighted graphs.
         Therefore all entries must be 0 or 1.")
  }

  if (! isSymmetric(adj_mat_1) | ! isSymmetric(adj_mat_2)){
    stop("Only undirected graphs are supported so arguments must be
         symmetric. This includes that rownames must equal colnames.")
  }

  if (any(c(diag(adj_mat_1), diag(adj_mat_2)) != 0)){
    stop("Loops are not supported so diagonal must be all 0.")
  }

  if (is.null(colnames(adj_mat_1)) | is.null(rownames(adj_mat_1))
      | is.null(colnames(adj_mat_2)) | is.null(rownames(adj_mat_2))){
    stop("The matrices must be named.")
  }

  if (! compare::compare(rownames(adj_mat_1), rownames(adj_mat_2),
                         ignoreOrder = TRUE)$result[1]){
    stop("The node names must be the same in both graphs.")
  }
  
  if (sum(adj_mat_1) != sum(adj_mat_2)){
    stop(paste("The graphs do not have the same number of edges indicating",
               "that the t-cherry trees do not have the same order", sep = " "))
  }
  
  adj_mat_1 <- adj_mat_1[rownames(adj_mat_2), colnames(adj_mat_2)]

  u_tri_1_idx <- upper.tri(adj_mat_1)
  u_tri_2_idx <- upper.tri(adj_mat_2)

  u_tri_1 <- adj_mat_1[u_tri_1_idx]
  u_tri_2 <- adj_mat_2[u_tri_2_idx]

  is_different <- u_tri_1 != u_tri_2
  length(which(is_different)) / 2
}
