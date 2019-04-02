#' Conditional independence test
#'
#' @description Makes a conditional independece test based on a
#' likelihood ratio test.
#'
#' @param var1,var2 Variables to test independence for.
#' @param cond Variables to condition on. If empty the test is just an
#' independence test.
#' @param data Data with realisations of the variables.
#' @param smooth Additional cell counts when estimating probabilities.
#' May be used to avoid zero probabilities.
#'
#' @details
#' The calculated test statistic is \deqn{
#' 2 \sum n_comb * log((P(var1, var2|cond)) /
#' (P(var1|cond)P(var2|cond)))
#' }
#' where the sum is over all states of the variables, n_comb is the
#' number of times the current combination of states has been seen in
#' the data and the probability distributions are estimated. This
#' estimation is done with maximum likelihood if \code{smooth} is set
#' to 0. If this gives zero probabilities the \code{smooth} argument
#' must be used to prevent it in order to calculate the test size.
#'
#' The test statistic follows assymptotically a chi-square distribution
#' with \eqn{(n1 - 1)(n2 - 1)n_cond} degrees of freedom, where \eqn{n1}
#' and \eqn{n2}
#' are the number of states for \code{var1} and \code{var2} and
#' \eqn{n_cond} is the number of configurations of the states of the
#' variables to be conditioned on.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item\code{chi_sq_statistic} The calculated test statistic.
#' \item\code{df} Degrees of freedom.
#' \item\code{p_value} The p-value of the test.
#' }
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
#' var1 <- c(sample(c(1, 2), 100, replace = TRUE))
#' var2 <- var1 + c(sample(c(1, 2), 100, replace = TRUE))
#' var3 <- var1 + c(sample(c(0, 1), 100, replace = TRUE,
#'                         prob = c(0.9, 0.1)))
#' var4 <- c(sample(c(1, 2), 100, replace = TRUE))
#' var5 <- var2 + var3
#'
#' data <- data.frame("var1" = as.character(var1),
#'                    "var2" = as.character(var2),
#'                    "var3" = as.character(var3),
#'                    "var4" = as.character(var4),
#'                    "var5" = as.character(var5))
#'
#'
#' cond_independence_test("var1", "var4", data = data, smooth = 0.1)
#' cond_independence_test("var2", "var3", cond = c("var1", "var5"),
#'                        data = data, smooth = 0.1)
#'
#' @export

cond_independence_test <- function(var1, var2, cond = c(), data,
                                   smooth = 0){

  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix.")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors.")
  }

  if (length(setdiff(c(var1, var2, cond), names(data))) != 0){
    stop("var1, var2 and the variables in cond must be variable names in data.")
  }

  if (length(smooth) > 1){
    stop("smooth must be a single non-negative value.")
  }
  else if (!is.numeric(smooth)) {
    stop("smooth must be numeric.")
  }
  else if (smooth < 0){
    stop("smooth must be a non-negative numeric value.")
  }

  if (length(cond) == 0){
    ncond <- 1

    tab_12 <- table(data[, c(var1, var2)])
    p_tab_12 <- (tab_12 + smooth) / sum(tab_12 + smooth)

    tab_1 <- table(data[, var1])
    names(dimnames(tab_1)) <- var1
    p_tab_1 <- (tab_1 + smooth) / sum(tab_1 + smooth)

    tab_2 <- table(data[, var2])
    names(dimnames(tab_2)) <- var2
    p_tab_2 <- (tab_2 + smooth) / sum(tab_2 + smooth)

    if (0 %in% c(p_tab_12, p_tab_1, p_tab_2)){
      stop("Some probabilities are zero. Consider using the smooth argument.")
    }

    chi_sq <- 2 *
      sum(gRbase::ar_prod(tab_12,
                          log(gRbase::ar_div(p_tab_12,
                                             gRbase::ar_prod(p_tab_1,
                                                             p_tab_2)))))
  }else{
    tab_cond <- table(data[, cond], dnn = cond)
    p_tab_cond <- (tab_cond + smooth) / sum(tab_cond + smooth)

    if (0 %in% c(p_tab_cond)){
      stop("Some probabilities are zero. Consider using the smooth argument.")
    }

    ncond <- length(tab_cond)

    tab_12cond <- table(data[, c(var1, var2, cond)])
    p_tab_12cond <- (tab_12cond + smooth) / sum(tab_12cond + smooth)
    p_tab_12_given_cond <- gRbase::ar_div(p_tab_12cond, p_tab_cond)

    tab_1cond <- table(data[, c(var1, cond)])
    p_tab_1cond <- (tab_1cond + smooth) / sum(tab_1cond + smooth)
    p_tab_1_given_cond <- gRbase::ar_div(p_tab_1cond, p_tab_cond)

    tab_2cond <- table(data[, c(var2, cond)])
    p_tab_2cond <- (tab_2cond + smooth) / sum(tab_2cond + smooth)
    p_tab_2_given_cond <- gRbase::ar_div(p_tab_2cond, p_tab_cond)

    if (0 %in% c(p_tab_12cond, p_tab_1cond, p_tab_2cond)){
      stop("Some probabilities are zero. Consider using the smooth argument.")
    }

    chi_sq <- 2 *
      sum(gRbase::ar_prod(tab_12cond,
                          log(gRbase::ar_div(
                            p_tab_12_given_cond,
                            gRbase::ar_prod(p_tab_1_given_cond,
                                            p_tab_2_given_cond)))))
  }

  n1 <- length(unique(data[, var1]))
  n2 <- length(unique(data[, var2]))

  df <- (n1 - 1) * (n2 - 1) * ncond

  p <- 1 - stats::pchisq(chi_sq, df)

  return(list("chi_sq_statistic" = chi_sq,
              "df" = df,
              "p_value" = p))

}
