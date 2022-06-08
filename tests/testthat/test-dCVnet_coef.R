
context("tests of coefficient extraction")

# see make_testdata.R for data/model.RData

load(file = "data/model.RData")
# m

# there are 30 predictors plus an intercept
#   evaluated for 2x10-fold repeated CV
# Thus 31 x 20 =

types_of_coefs <- c("production",
                    "all",
                    "mean",
                    "median",
                    "byrep",
                    "byrep_mean",
                    "byrep_median")
names(types_of_coefs) <- types_of_coefs
coef_list <- lapply(types_of_coefs,
                    function(x) coef(m, type = x))

coef_nrow <- vapply(coef_list, nrow, FUN.VALUE = 1L)
coef_ncol <- vapply(coef_list, ncol, FUN.VALUE = 1L)

expected_nrow <- 31L * c(1L, 20L, 1L, 1L, 2L, 1L, 1L)
expected_ncol <- c(2L, 4L, 2L, 2L, 3L, 2L, 2L)
names(expected_ncol) <- names(expected_nrow) <- types_of_coefs

# Use a precalculated model (see above):
test_that(
  "Coef produces the expected dimensions", {
  # check rows:
  expect_equal(coef_nrow, expected_nrow)
  # check cols:
  expect_equal(coef_ncol, expected_ncol)
  }
)

test_that(
  "nonsensical coef types produce errors",
  expect_error(coef(m, type = "hafdsa"))
)
