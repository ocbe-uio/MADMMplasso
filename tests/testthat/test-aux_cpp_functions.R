x <- 1:12
test_that("multiples_of() works", {
  expect_equal(multiples_of(x, 1), matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
  expect_equal(multiples_of(x, 2), matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)))
  expect_equal(multiples_of(x, 3), matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1)))
  expect_equal(multiples_of(x, 4), matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)))
  expect_equal(multiples_of(x, 5), matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0)))
  expect_equal(multiples_of(x, 6), matrix(c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1)))
})
