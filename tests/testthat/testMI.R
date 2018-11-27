context('MI2')
library(tcherry)

test_that('Input is correctly specified',{
  expect_error(MI2(c(1.5,3.6),c('1','2')),'x and y must be either characters or factors')
  expect_error(MI2(c('1','2'),c('1','2'),c(1,2)),'smooth must be a single positive number')
  expect_error(MI2(c('1','2'),c('1','2'),'C'),'smooth must be numeric')
  expect_error(MI2(c('1','2'),c('1','2'),-1),'smooth must be a non-negative numeric value')
})


x <- c('hund','kat','kat','kat','hund','hund','hund','hund','hund')
y <- c('a','b','c','a','c','b','a','a','b')

tab_x <- table(x)
tab_y <- table(y)
tab_xy <- table(x,y)

prop_x <- c(tab_x/sum(tab_x))
prop_y <- c(tab_y/sum(tab_y))
prop_xy <- tab_xy/sum(tab_xy)

frac_prop_MI <- t(t(prop_xy)/prop_y)/prop_x

result <- as.table(matrix(c(1.125, 0.750, 1, 1, 0.75, 1.5),nr=2))
dimnames(result) <- list('x'=c('hund','kat'),'y'=c('a','b','c'))

test_that('frac for MI is calculated properly',{
  expect_equal(frac_prop_MI,result)
})

log_frac <- log(frac_prop_MI, 2)

log_frac_result <- as.table(matrix(c(0.1699250, 0, -0.4150375, -0.4150375, 0, 0.5849625),nr=2, byrow = TRUE))
dimnames(log_frac_result) <- list('x'=c('hund','kat'),'y'=c('a','b','c'))

prod_prop <- prop_xy*log_frac

prod_prop_result <- as.table(matrix(c(0.05664167, 0, -0.04611528, -0.04611528, 0, 0.06499583),nr=2, byrow = TRUE))
dimnames(prod_prop_result) <- list('x'=c('hund','kat'),'y'=c('a','b','c'))

MI_result <- 0.02940695

test_that('MI is correct',{
  expect_equal(log_frac, log_frac_result)
  expect_equal(prod_prop, prod_prop_result, tolerance = 10^{-5})
  expect_equal(MI2(x,y),MI_result)
})

x <- c('1','2','3')
y <- c('a','b','b')

test_that('an error occur when zero probabilities',{
  expect_error(MI2(x,y), 'Some probabilities are zero and therefore MI cannot be calculated. Consider using the smooth argument.')
})


