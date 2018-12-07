#' Estimate conditional probability tables
#'
#' @description Estimates the conditional probability tables for
#' bayesian network models, where the structure is given by an
#' adjacency matrix.
#'
#' @param adj_matrix The adjacency matrix for the DAG.
#' @param data The data the probabilities should be estimated from.
#' @param bayes_smooth The additional cell counts for
#' bayesian estimation.
#'
#' @return A list of the conditional probability tables for the
#' bayesian network. If the \code{bayes_smooth} argument is zero,
#' it is the maximum likelihood estimates. Otherwise, it is bayesian
#' estimates.
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @examples
#' set.seed(43)
#' var1 <- c(sample(c(1, 2), 50, replace = TRUE))
#' var2 <- var1 + c(sample(c(1, 2), 50, replace = TRUE))
#' var3 <- var1 + c(sample(c(0, 1), 50, replace = TRUE,
#'                         prob = c(0.9, 0.1)))
#' var4 <- c(sample(c(1, 2), 50, replace = TRUE))
#' data <- data.frame("var1" = as.character(var1),
#'                    "var2" = as.character(var2),
#'                    "var3" = as.character(var3),
#'                    "var4" = as.character(var4))
#'
#' adj_matrix_DAG <- matrix(c(0, 0, 0, 0,
#'                            1, 0, 0, 0,
#'                            1, 0, 0, 0,
#'                            0, 1, 0, 0),
#'                           nrow = 4)
#' CPT(adj_matrix_DAG, data)
#' CPT(adj_matrix_DAG, data, bayes_smooth = 1)
#' @export

CPT <- function(adj_matrix, data, bayes_smooth = 0){
  data <- data.frame(data, stringsAsFactors = FALSE)
  nodes <- rownames(adj_matrix)
  FUN <- function(node){
    parents_idx <- which(adj_matrix[, node] == 1)
    parents <- nodes[parents_idx]

    tab <- table(data[, c(node, parents)]) + bayes_smooth
    if (length(parents) == 0){
      mar <- NULL
      names(dimnames(tab)) <- node
    } else {
      tab_parents <- table(data[, c(parents)]) + bayes_smooth
      if (any(tab_parents == 0)){
        stop("Some cell counts of parent configurations are zero.
             Consider using the bayes_smooth argument.")
      }
      mar <- (1:length(parents)) + 1
      }
    prop.table(tab, margin = mar)
  }

  CPT_list <- lapply(nodes, FUN)
  names(CPT_list) <- nodes
  CPT_list
}
