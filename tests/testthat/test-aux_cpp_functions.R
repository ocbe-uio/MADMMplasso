x <- 1:12
set.seed(1282897)
y <- rpois(n = 10, lambda = 10)
test_that("multiples_of() works", {
  expect_error(multiples_of(x, 0), "divisor cannot be 0")
  expect_identical(multiples_of(x, 1), matrix(1:12))
  expect_identical(multiples_of(x, 2), matrix(c(2L, 4L, 6L, 8L, 10L, 12L)))
  expect_identical(multiples_of(x, 3), matrix(c(3L, 6L, 9L, 12L)))
  expect_identical(multiples_of(x, 4), matrix(c(4L, 8L, 12L)))
  expect_identical(multiples_of(x, 5), matrix(c(5L, 10L)))
  expect_identical(multiples_of(x, 6), matrix(c(6L, 12L)))

  expect_identical(multiples_of(y, 4), matrix(c(2L, 3L, 4L, 6L)))
  expect_identical(multiples_of(y, 5), matrix(c(1L, 5L, 9L)))
  expect_identical(multiples_of(y, 7), matrix(c(8L, 10L)))
  expect_identical(multiples_of(y, 15), matrix(c(1L, 5L, 9L)))
  expect_identical(multiples_of(y, 10), matrix(as.integer(NULL)))

  expect_identical(multiples_of(y, 4, TRUE), matrix(c(y[2], y[3], y[4], y[6])))
  expect_identical(multiples_of(y, 5, TRUE), matrix(c(y[1], y[5], y[9])))
  expect_identical(multiples_of(y, 7, TRUE), matrix(c(y[8], y[10])))
  expect_identical(multiples_of(y, 15, TRUE), matrix(c(y[1], y[5], y[9])))
  expect_identical(multiples_of(y, 10, TRUE), matrix(as.integer(NULL)))
})

N <- rpois(1, 10)
p <- rpois(1, 10)
x <- matrix(rnorm(N * p), N, p)
y <- matrix(rnorm(p), p)

test_that("C++ scale() works like R scale()", {
  expect_equal(scale(x, FALSE, y), scale_cpp(x, y), ignore_attr = TRUE)
  expect_equal(scale(x, FALSE, 1 / y), scale_cpp(x, 1 / y), ignore_attr = TRUE)
})

test_that("sqrt_sum_squared_rows() works", {
  sssr <- function(x) sqrt(apply(x^2, 1, sum, na.rm = TRUE))
  expect_equal(sssr(x), sqrt_sum_squared_rows(x), ignore_attr = TRUE)
  expect_equal(sssr(y), sqrt_sum_squared_rows(y), ignore_attr = TRUE)
  expect_equal(sssr(t(y)), sqrt_sum_squared_rows(t(y)), ignore_attr = TRUE)
})

test_that("modulo() works", {
  top <- rpois(1, 20)
  for (i in seq_len(top)) {
    x <- matrix(rpois(top, i))
    expect_equal(x %% i, modulo(x, i), ignore_attr = TRUE)
  }
})

test_that("count_nonzero_a() works", {
  dims <- rpois(5L, 10L)
  x <- matrix(rpois(prod(dims[1:2]), 1), dims[1], dims[2])
  y <- matrix(rpois(prod(dims[1:2]), 1), dims[1], dims[2])
  z <- array(rpois(prod(dims[3:5]), 1), dims[3:5])
  w <- array(rpois(prod(dims[3:5]), 1), dims[3:5])
  expect_equal(count_nonzero_a(x), count_nonzero_a_cpp(x), ignore_attr = TRUE)
  expect_equal(count_nonzero_a(y), count_nonzero_a_cpp(y), ignore_attr = TRUE)
  expect_equal(count_nonzero_a(z), count_nonzero_a_cpp(z), ignore_attr = TRUE)
  expect_equal(count_nonzero_a(w), count_nonzero_a_cpp(w), ignore_attr = TRUE)
  expect_equal(count_nonzero_a(x), count_nonzero_a_sp_mat(as(x, "sparseMatrix")), ignore_attr = TRUE)
  expect_equal(count_nonzero_a(y), count_nonzero_a_sp_mat(as(y, "sparseMatrix")), ignore_attr = TRUE)
  expect_equal(count_nonzero_a(z), count_nonzero_a_cube(z), ignore_attr = TRUE)
  expect_equal(count_nonzero_a(w), count_nonzero_a_cube(w), ignore_attr = TRUE)
})
