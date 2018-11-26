context('sum')
library(tcherry)

a <- 7
b <- 5
test_that('addition is right',{
  expect_equal(sum(1,2),3)
  expect_equal(sum(2,2),4)
})
