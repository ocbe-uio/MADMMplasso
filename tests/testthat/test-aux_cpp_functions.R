x <- 1:12
set.seed(1282897)
y <- rpois(n = 10, lambda = 10)
test_that("multiples_of() works", {
  expect_error(multiples_of(x, 0), "divisor cannot be 0")
  expect_equal(multiples_of(x, 1), matrix(1:12))
  expect_equal(multiples_of(x, 2), matrix(c(2, 4, 6, 8, 10, 12)))
  expect_equal(multiples_of(x, 3), matrix(c(3, 6, 9, 12)))
  expect_equal(multiples_of(x, 4), matrix(c(4, 8, 12)))
  expect_equal(multiples_of(x, 5), matrix(c(5, 10)))
  expect_equal(multiples_of(x, 6), matrix(c(6, 12)))

  expect_equal(multiples_of(y, 4), matrix(c(2, 3, 4, 6)))
  expect_equal(multiples_of(y, 5), matrix(c(1, 5, 9)))
  expect_equal(multiples_of(y, 7), matrix(c(8, 10)))
  expect_equal(multiples_of(y, 15), matrix(c(1, 5, 9)))
  expect_equal(multiples_of(y, 10), matrix(as.integer(NULL)))

  expect_equal(multiples_of(y, 4, TRUE), matrix(c(y[2], y[3], y[4], y[6])))
  expect_equal(multiples_of(y, 5, TRUE), matrix(c(y[1], y[5], y[9])))
  expect_equal(multiples_of(y, 7, TRUE), matrix(c(y[8], y[10])))
  expect_equal(multiples_of(y, 15, TRUE), matrix(c(y[1], y[5], y[9])))
  expect_equal(multiples_of(y, 10, TRUE), matrix(as.integer(NULL)))
})
