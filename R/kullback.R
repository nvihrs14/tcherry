#' Calculates the Kullback-Leibler divergence
#'
#' @description Calculates the Kullback-Leibler divergence for two
#' discrete probability distributions over the same universe.
#'
#' @param target_p Named array with the target distribution.
#' @param approx_p Named array with the approximate distribution.
#' @param log_base The base used for the logarithm.
#'
#' @details The Kullback-Leibler divergence between two probability
#' distributions p and q is calculated by the
#' formula \deqn{D(p||q) = \sum p(x)log(p(x)/q(x))} where
#' the sum is over all possible states of x.
#'
#' @return The Kullback-Leibler divergence.
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @references
#' \insertRef{TCJT}{tcherry}
#'
#' @examples
#' p_target <- array(c(0.1, 0.2, 0.05, 0.4, 0.005, 0.015, 0.13, 0.1),
#'                   dim = c(2, 2, 2),
#'                   dimnames = list("V1" = c("a", "b"),
#'                                   "V2" = c("a", "b"),
#'                                   "V3" = c("a", "b")))
#'                                   
#' p_approx <- array(c(0.01, 0.23, 0.06, 0.4, 0.005, 0.018,
#'                     0.137, 0.14),
#'                   dim = c(2, 2, 2),
#'                   dimnames = list("V1" = c("a", "b"),
#'                                   "V2" = c("a", "b"),
#'                                   "V3" = c("a", "b")))
#'                                   
#' kullback(p_target, p_approx)
#' @export

kullback <- function(target_p, approx_p, log_base = 2){
  
  if (! is.array(target_p) | !is.array(approx_p)){
    stop("The input is not arrays.")
  }

  if (is.null(dimnames(target_p)) | is.null(dimnames(approx_p))){
    stop("The arrays must be named.")
  }

  if (is.null(names(dimnames(target_p))) |
      is.null(names(dimnames(approx_p)))){
    stop("Variable names are missing from the arrays.")
  }

  if (0 %in% lengths(dimnames(target_p)) |
      0 %in% lengths(dimnames(approx_p))){
    stop("Names of the variable levels are missing.")
  }

  if (!isTRUE(all.equal(c(sum(target_p), sum(approx_p)), c(1, 1)))){
    stop("At least one of the given arrays is not a probability distribution.")
  }

  if (0 %in% target_p | 0 %in% approx_p){
    stop("Some probabilities are zero.")
  }

  if (! compare::compare(dimnames(target_p), dimnames(approx_p),
                       ignoreComponentOrder = TRUE)$result[1]){
    stop("Distributions are not over the same universe.")
  }
  
  frac <- gRbase::ar_div(target_p, approx_p)
  sum(gRbase::ar_prod(target_p, log(frac, base = log_base)))
}