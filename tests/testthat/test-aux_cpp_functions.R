x <- 1:12
test_that("multiples_of() works", {
  expect_error(multiples_of(x, 0), "divisor cannot be 0")
  expect_equal(multiples_of(x, 1), matrix(1:12))
  expect_equal(multiples_of(x, 2), matrix(c(2, 4, 6, 8, 10, 12)))
  expect_equal(multiples_of(x, 3), matrix(c(3, 6, 9, 12)))
  expect_equal(multiples_of(x, 4), matrix(c(4, 8, 12)))
  expect_equal(multiples_of(x, 5), matrix(c(5, 10)))
  expect_equal(multiples_of(x, 6), matrix(c(6, 12)))
})

# TODO: test indices, instead of sequence
