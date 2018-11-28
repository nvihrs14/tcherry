

#' Calculate mutual information
#' @description Calculate mutual information for two or three            categorical variables.
#'
#' @param x,y,z Vectors of class character or factor.
#' @param smooth A single numeric non-negative value. See "details"      below.
#' @param log_base The base of the logarithmic function to be used.
#' @details The \code{smooth} argument is added to all cell counts in     tables before normalisation. This is meant as a way to avoid zero     counts, leading to zero probabilities.
#'
#' The mutual information for two variables is calculated by the         formula \deqn{MI(x, y) = \sum P(x, y) log(P(x, y) / (P(x)P(y)))}     where the sum is over alle possible values of x and y.
#'
#' The mutual information for three variables is calculated by the      formula \deqn{MI(x, y, z) = \sum P(x, y, z) log(P(x, y, z) / (P      (x)P(y)P(z)))} where the sum is over all possible values of x, y     and z.
#' @return The mutual information given by a single numeric value.
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{TCJT}{tcherry}
#'
#' \insertRef{EKTS}{tcherry}
#' @examples
#'
#' x <- c("hund","kat","kat","kat","hund","hund","hund","hund","hund")
#' y <- c("a","b","c","a","c","b","a","a","b")
#' MI2(x, y)
#' MI2(x, y, smooth = 1)
#' MI2(x, y, log_base = exp(1))
#' @export
MI2 <- function(x, y, smooth = 0, log_base = 2){
  if (! all(sapply(list(x, y), function(x){
    is.character(x) | is.factor(x)
    }
    ))){
    stop("x and y must be either characters or factors")
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

  tab_x <- table(x) + smooth
  tab_y <- table(y) + smooth
  tab_xy <- table(x, y) + smooth

  if (0 %in% c(tab_x, tab_xy, tab_y)){
    stop("Some probabilities are zero and therefore MI cannot be calculated.
         Consider using the smooth argument.")
  }

  prop_x <- c(tab_x / sum(tab_x))
  prop_y <- c(tab_y / sum(tab_y))
  prop_xy <- tab_xy / sum(tab_xy)

  frac_prop_MI <- sweep(prop_xy, 1, prop_x, FUN = "/")
  frac_prop_MI <- sweep(frac_prop_MI, 2, prop_y, FUN = "/")

  MI <- sum(prop_xy * log(frac_prop_MI, base = log_base))

  return(MI)
}

#' @rdname MI2
MI3 <- function(x, y, z, smooth = 0, log_base = 2){
  if (! all(sapply(list(x, y, z), function(x){
    is.character(x) | is.factor(x)
    }
    ))){
    stop("x, y and z must be either characters or factors")
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

  tab_x <- table(x) + smooth
  tab_y <- table(y) + smooth
  tab_z <- table(z) + smooth
  tab_xyz <- table(x, y, z) + smooth

  if (0 %in% c(tab_x, tab_y, tab_z, tab_xyz)){
    stop("Some probabilities are zero and therefore MI cannot be calculated.
         Consider using the smooth argument.")
  }

  prop_x <- c(tab_x / sum(tab_x))
  prop_y <- c(tab_y / sum(tab_y))
  prop_z <- c(tab_z / sum(tab_z))
  prop_xyz <- tab_xyz / sum(tab_xyz)

  frac_prop_MI <- sweep(prop_xyz, 1, prop_x, FUN = "/")
  frac_prop_MI <- sweep(frac_prop_MI, 2, prop_y, FUN = "/")
  frac_prop_MI <- sweep(frac_prop_MI, 3, prop_z, FUN = "/")

  MI <- sum(prop_xyz * log(frac_prop_MI, base = log_base))

  return(MI)
}
