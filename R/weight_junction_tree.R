#' Calculates the weight of a junction tree
#' 
#' @description Calculates the weight of a junction tree with the given
#' cliques and separators.
#'
#' @param cliques,separators Lists containing cliques and separators
#' given by variable names.
#' @param data A data frame or matrix containing the variables.
#' @param ... Additional arguments passed to \code{MIk}.
#'
#' @details The weight of the junction tree is calculated
#' as \deqn{\sum MI(clique)-\sum MI(separator)} where the sum is over
#' all cliques and separators respectively.
#' 
#' @return The weight of the junction tree.
#' 
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#' 
#' @importFrom Rdpack reprompt
#' 
#' @references
#' \insertRef{EKTShyp}{tcherry}
#' 
#' @seealso
#' \code{\link{MIk}} for mutual infomation for k variables.
#' 
#' @examples
#' set.seed(43)
#' var1 <- c(sample(c(1, 2), 100, replace = TRUE))
#' var2 <- var1 + c(sample(c(1, 2), 100, replace = TRUE))
#' var3 <- var1 + c(sample(c(0, 1), 100, replace = TRUE,
#'                         prob = c(0.9, 0.1)))
#' var4 <- c(sample(c(1, 2), 100, replace = TRUE))
#'
#' data <- data.frame("var1" = as.character(var1),
#'                    "var2" = as.character(var2),
#'                    "var3" = as.character(var3),
#'                    "var4" = as.character(var4))
#'                    
#' cliques <- list(c("var1", "var2", "var4"), c("var1", "var3", "var4"))
#' separators <- list(c("var1", "var4"))
#'
#' weight_junction_tree(cliques, separators, data, smooth = 0.001)
#' @export

weight_junction_tree <- function(cliques, separators, data, ...){

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

  if (!all(c(unlist(cliques), unlist(separators)) %in% colnames(data))){
    stop("All names in cliques and separators must be column names of data.")
  }

  if (length(cliques) - 1 != length(separators)){
    stop("The number of separators must be one less than the number of cliques.
         If a separator is used more than once it should be repeated in the list.")
  }

  MI_cliq <- lapply(cliques, MIk, data, ...)
  MI_sep <- lapply(separators, MIk, data, ...)

  sum(unlist(MI_cliq)) - sum(unlist(MI_sep))
}