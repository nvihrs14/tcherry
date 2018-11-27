#' Calculate mutual information
#' @description Calculate mutual information for two or three variables.
#'
#' @param x,y,z Vectors of class character or factor.
#' @param smooth A single numeric non-negative value. See 'details' below.
#' @param log_base The base of the logarithmic function to be used.
#' @details The \code{smooth} argument is added to all cell counts in tables before normalisation. This is meant as a way to avoid zero counts, leading to zero probabilities.
#'
#'     The mutual information for two variables is calculated by the formula \deqn{MI(x,y)=\sum P(x,y)log(P(x,y)/P(x)P(y))} where the sum is over alle possible values of x and y.
#' @return The mutual information given by a single numeric value.
MI2 <- function(x, y, smooth = 0, log_base = 2){
  if (! all(sapply(list(x, y), function(x){is.character(x) | is.factor(x)}))){
    stop('x and y must be either characters or factors')
  }
  if (length(smooth) > 1){
    stop('smooth must be a single non-negative value')
  }
  else if (!is.numeric(smooth)) {
    stop('smooth must be numeric')
  }
  else if (smooth < 0){
    stop('smooth must be a non-negative numeric value')
  }

  tab_x <- table(x) + smooth
  tab_y <- table(y) + smooth
  tab_xy <- table(x, y) + smooth

  if (0 %in% c(tab_x, tab_xy, tab_y)){
    stop('Some probabilities are zero and therefore MI cannot be calculated. Consider using the smooth argument.')
  }

  prop_x <- c(tab_x / sum(tab_x))
  prop_y <- c(tab_y / sum(tab_y))
  prop_xy <- tab_xy / sum(tab_xy)

  frac_prop_MI <- t(t(prop_xy) / prop_y) / prop_x

  MI <- sum(prop_xy * log(frac_prop_MI, base = log_base))

  return(MI)
}


