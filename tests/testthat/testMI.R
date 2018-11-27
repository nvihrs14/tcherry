context('MI2')
library(tcherry)

test_that('Input is correctly specified',{
  expect_error(MI2(c(1.5,3.6),c('1','2')),'x and y must be either characters or factors')
  expect_error(MI2(c('1','2'),c('1','2'),c(1,2)),'smooth must be a single non-negative value')
  expect_error(MI2(c('1','2'),c('1','2'),'C'),'smooth must be numeric')
  expect_error(MI2(c('1','2'),c('1','2'),-1),'smooth must be a non-negative numeric value')
})


x <- c('hund','kat','kat','kat','hund','hund','hund','hund','hund')
y <- c('a','b','c','a','c','b','a','a','b')

MI_result <- 0.02940695

test_that('MI is correct',{
  expect_equal(MI2(x,y),MI_result)
})

x1 <- c('1','2','3')
y1 <- c('a','b','b')

test_that('an error occur when zero probabilities', {
  expect_error(MI2(x1, y1), 'Some probabilities are zero and therefore MI cannot be calculated. Consider using the smooth argument.')
})

context('MI3')

test_that('Input is correctly specified',{
  expect_error(MI3(c(1.5, 3.6), c('1', '2'), c('a', 'b')),'x, y and z must be either characters or factors')
  expect_error(MI3(c('1', '2'), c('1', '2'), c('1', '2'), c(1, 2)),'smooth must be a single non-negative value')
  expect_error(MI3(c('1', '2'), c('1', '2'), c('1', '2'), 'C'),'smooth must be numeric')
  expect_error(MI3(c('1', '2'), c('1', '2'), c('1', '2'), -1),'smooth must be a non-negative numeric value')
})

x1 <- c('1', '2', '3')
y1 <- c('a', 'b', 'b')
z1 <- c('hej', 'goddag', 'hej')

test_that('an error occur when zero probabilities', {
  expect_error(MI3(x1, y1, z1), 'Some probabilities are zero and therefore MI cannot be calculated. Consider using the smooth argument.')
})

x <- c('hund','kat','kat','kat','hund','hund','hund','hund','hund')
y <- c('a','b','c','a','c','b','a','a','b')
z <- c('1', '2', '1', '1', '2', '1', '2', '2', '2')
smooth <- 0.1

MI_result_3 <- 0.3172042

test_that('MI is calculated correctly', {
  expect_equal(MI3(x, y, z, smooth = smooth), MI_result_3)
})
