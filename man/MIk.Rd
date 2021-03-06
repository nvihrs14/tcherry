% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MI_k.R
\name{MIk}
\alias{MIk}
\title{Calculates mutual information}
\usage{
MIk(variables, data, smooth = 0, log_base = 2)
}
\arguments{
\item{variables}{Vector of variable names.}

\item{data}{Data frame or matrix with the variables.}

\item{smooth}{Additional cell counts for bayesian estimation of
probabilities.}

\item{log_base}{The base of the logarithmic function to be used.}
}
\value{
The mutual information given by a single numeric value.
}
\description{
Calculates mutual information for k categorical
variables.
}
\details{
The mutual information for variables \eqn{V1,\ldots,Vk} is calculated
by the formula \deqn{MI(V1,\ldots,Vk) = \sum P(V1,\ldots,Vk) log(P
(V1,\ldots,Vk) / (P(V1)\dots P(Vk)))} where the sum is over all
possible values of \eqn{V1,\ldots,Vk}.

If the function is used for two variables it corresponds to using
\code{MI2} and for three variables it corresponds to \code{MI3}.
}
\examples{
set.seed(43)
var1 <- c(sample(c(1, 2), 100, replace = TRUE))
var2 <- var1 + c(sample(c(1, 2), 100, replace = TRUE))
var3 <- var1 + c(sample(c(0, 1), 100, replace = TRUE,
                        prob = c(0.9, 0.1)))
var4 <- c(sample(c(1, 2), 100, replace = TRUE))
var5 <- var2 + var3
var6 <- var1 - var4 + c(sample(c(1, 2), 100, replace = TRUE))
var7 <- c(sample(c(1, 2), 100, replace = TRUE))

data <- data.frame("var1" = as.character(var1),
                   "var2" = as.character(var2),
                   "var3" = as.character(var3),
                   "var4" = as.character(var4),
                   "var5" = as.character(var5),
                   "var6" = as.character(var6),
                   "var7" = as.character(var7))
                   
MIk(c("var1", "var2", "var7"), data, smooth = 0.001)
}
\references{
\insertRef{EKTShyp}{tcherry}
}
\seealso{
\code{\link{MI2}} and \code{\link{MI3}} for mutual infomation for two
or three variables respectively.
}
\author{
Katrine Kirkeby, \email{enir_tak@hotmail.com}

Maria Knudsen, \email{mariaknudsen@hotmail.dk}

Ninna Vihrs, \email{ninnavihrs@hotmail.dk}
}
