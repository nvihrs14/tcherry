# Hello, world!  This is an example function named 'hello' which prints
# 'Hello, world!'.  You can learn more about package authoring with RStudio
# at: http://r-pkgs.had.co.nz/ Some useful keyboard shortcuts for package
# authoring: Build and Reload Package: 'Ctrl + Shift + B' Check Package:
# 'Ctrl + Shift + E' Test Package: 'Ctrl + Shift + T'

hello <- function(){
  print("Hello, world!")
}

#' Add two numbers
#'
#' @param x,y numeric numbers
#' @return The sum of \code{x} and \code{y}
#' @examples
#' sum(2,3)
sum <- function(x,y){
  x+y
}

