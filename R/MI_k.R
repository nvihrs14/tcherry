#' Calculate mutual information
#' @description Calculate mutual information for k categorical variables.
#'
#' @param variables Vector of variable names.
#' @param data Data frame or matrix with the variables.
#' @param smooth A single numeric non-negative value. See "details"      below.
#' @param log_base The base of the logarithmic function to be used.
#' @details The \code{smooth} argument is added to all cell counts in     tables before normalisation. This is meant as a way to avoid zero     counts, leading to zero probabilities.
#'
#' The mutual information for variables \eqn{V1,\ldots,Vk} is calculated
#' by the formula \deqn{MI(V1,\ldots,Vk) = \sum P(V1,\ldots,Vk) log(P(V1,\ldots,Vk) / (P(V1)\dots P(Vk)))} where the sum is over all
#' possible values of \eqn{V1,\ldots,Vk}.
#'
#' If the function is used for two variables it corresponds to using
#' \code{MI2} and for three variables it corresponds to \code{MI3}.
#' @return The mutual information given by a single numeric value.
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{EKTShyp}{tcherry}
#' @seealso
#' \code{\link{MI2}} and \code{\link{MI3}} for mutual infomation for two
#' or three variables respectively.
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
#' MIk(c("var1", "var2", "var7"), data, smooth = 0.001)
#' @export
MIk <- function(variables, data, smooth = 0, log_base = 2){
  if (!all(variables %in% colnames(data))){
    stop("All names in variables must be column names of data")
  }

  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix")
  }

  if (! all(sapply(as.list(data[, variables]), function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Input must be either characters or factors")
  }

  if (length(smooth) > 1){
    stop("smooth must be a single non-negative value")
  }
  else if (!is.numeric(smooth)) {
    stop("smooth must be numeric")
  }
  else if (smooth < 0){
    stop("smooth must be a non-negative numeric value")
  }

  tab_marginals <- lapply(as.list(variables), function(x){
    t <- table(data[, x]) + smooth
    names(dimnames(t)) <- x
    t
    }
    )
  tab_joint <- table(data[, variables]) + smooth

  if (0 %in% c(unlist(tab_marginals), tab_joint)){
    stop("Some probabilities are zero and therefore MI cannot be calculated.
         Consider using the smooth argument.")
  }

  prop_marginals <- lapply(tab_marginals, function(tab){
    tab / sum(tab)
  }
  )
  prop_joint <- tab_joint / sum(tab_joint)

  frac_prop_MI <- gRbase::ar_div(prop_joint,
                                 gRbase::ar_prod_list(prop_marginals))

  MI <- sum(gRbase::ar_prod(prop_joint,
                            log(frac_prop_MI, base = log_base)))

  return(MI)
  }
