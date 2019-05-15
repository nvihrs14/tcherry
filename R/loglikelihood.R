entropy <- function(var, smooth = 0, log_base = 2){
  tab <- table(var) + smooth
  p_tab <- tab / sum(tab)

  - sum(p_tab * log(p_tab, base = log_base))
}

#' Calculates the log-likelihood value
#'
#' @description Calculates the value of the log-likelihood function for a
#' graphical model with the given junction tree after inserting estimated
#' probabilities.
#'
#' @param cliques A list containing the cliques of the junction tree.
#' @param separators A list containing the separators of the junction
#' tree.
#' @param data A data frame with observations of the variables.
#' @param ... Additional arguments passed to \code{weight_junction_tree}.
#'
#' @details Calculates the value as \deqn{n_obs * (weight - \sum H(A))}
#'  where \eqn{n_obs} is the number of observations in data, \eqn{weight}
#'  is the estimated weight of the junction tree and \eqn{H} is the
#'  estimated entropy of the variable \eqn{A}. The sum is over all
#'  variables in data. All estimated probabilities are maximum likelihood
#'  estimates unless the smooth argument has been used to deal with zero
#'  probabilities.
#'
#' @return The calculated log-likelihood value.
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @seealso \code{\link{weight_junction_tree}} for calculating the weight
#' of a junction tree.
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
#' # smooth is used to deal with zero probabilities.
#' loglikelihood(cliques, separators, data, smooth = 0.1)
#' @export

loglikelihood <- function(cliques, separators, data, ...){

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

  n_obs <- nrow(data)

  weight <- weight_junction_tree(cliques, separators, data, ...)

  sum_entropy <- sum(sapply(data, entropy, ...))

  n_obs * (weight - sum_entropy)
}
