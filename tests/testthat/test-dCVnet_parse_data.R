
nested_anyna <- function(x) {
  any(vapply(x, function(k) any(is.na(k)), FUN.VALUE = FALSE))
}

# init some data:
y1 <- rep(0:1, each = 10)
set.seed(42)
x1 <- data.frame(matrix(rnorm(20 * 10), nrow = 20, ncol = 10))

# Versions with missing data:
y2 <- y1
y2[c(8, 12)] <- NA

x2 <- x1
x2[sample(NROW(x2), size = 3), sample(NCOL(x2), size = 3)] <- NA


# minimal call:
p1.basic <- dCVnet::parse_dCVnet_input(y = y1, data = x1, family = "binomial")

# yname:
yn <- paste0(sample(letters, size = 7), collapse = "")
p1.yname <- dCVnet::parse_dCVnet_input(y = y1, data = x1, family = "binomial",
                                       yname = yn)


# basic tests -------------------------------------------------------------

test_that("parsed data equals input for simple calls", {
  # outcome values:
  expect_equal(as.factor(y1), p1.basic$y)
  expect_equal(as.factor(y1), p1.yname$y)
  # data values:
  expect_equal(as.matrix(x1), p1.basic$x_mat, ignore_attr = TRUE)
  expect_equal(as.matrix(x1), p1.yname$x_mat, ignore_attr = TRUE)
  # data colnames:
  expect_equal(colnames(p1.basic$x_mat), colnames(x1))
  expect_equal(colnames(p1.yname$x_mat), colnames(x1))

  # yname:
  expect_equal("y", p1.basic$yname)
  expect_equal(yn, p1.yname$yname)

  # output when yname changes:
  expect_equal(p1.basic[names(p1.basic) != "yname"],
               p1.yname[names(p1.yname) != "yname"])
})


# missing data ------------------------------------------------------------

test_that("parsed data works with missing data", {
  # First y only missing:
  expect_warning((p2.y <- dCVnet::parse_dCVnet_input(y = y2,
                                                     data = x1,
                                                     family = "binomial")))
  expect_equal(NROW(p2.y$y), NROW(p2.y$x_mat))
  expect_equal(NROW(p2.y$y), sum(complete.cases(y2)))
  expect_equal(nested_anyna(p2.y), FALSE)
  # Second x only missing:
  expect_warning((p2.x <- dCVnet::parse_dCVnet_input(y = y1,
                                                     data = x2,
                                                     family = "binomial")))
  expect_equal(NROW(p2.x$y), NROW(p2.x$x_mat))
  expect_equal(NROW(p2.x$y), sum(complete.cases(x2)))
  expect_equal(nested_anyna(p2.x), FALSE)

  # Third x&y missing:
  expect_warning((p2.xy <- dCVnet::parse_dCVnet_input(y = y2,
                                                      data = x2,
                                                      family = "binomial")))
  expect_equal(NROW(p2.xy$y), NROW(p2.xy$x_mat))
  expect_equal(NROW(p2.xy$y), sum(complete.cases(x2) & complete.cases(y2)))
  expect_equal(nested_anyna(p2.xy), FALSE)
})

# as above, but now passNA:
test_that("missing x can be passed (but not y)", {
  # First y only missing (no effect of passNA):
  expect_warning((p2.y <- dCVnet::parse_dCVnet_input(y = y2,
                                                     data = x1,
                                                     family = "binomial",
                                                     passNA = FALSE)))
  expect_warning((p3.y <- dCVnet::parse_dCVnet_input(y = y2,
                                                     data = x1,
                                                     family = "binomial",
                                                     passNA = TRUE)))

  expect_identical(p3.y, p2.y)
  # Second x only missing:
  (p3.x <- dCVnet::parse_dCVnet_input(y = y1,
                                      data = x2,
                                      family = "binomial",
                                      passNA = TRUE))
  expect_equal(NROW(p3.x$y), NROW(p3.x$x_mat))
  expect_gt(NROW(p3.x$y), sum(complete.cases(x2)))
  expect_equal(nested_anyna(p3.x), TRUE)

  # Third x&y missing:

  suppressWarnings(expect_warning((
    p3.xy <- dCVnet::parse_dCVnet_input(
      y = y2,
      data = x2,
      family = "binomial",
      passNA = TRUE
    )
  )))
  expect_equal(NROW(p3.xy$y), NROW(p3.xy$x_mat))
  expect_gt(NROW(p3.xy$y), sum(complete.cases(x2) & complete.cases(y2)))
  expect_equal(nested_anyna(p3.xy), TRUE)
})


# Family Handling ---------------------------------------------------------

# ~ binomial --------------------------------------------------------------

# make some factor labels
set.seed(42)
y_fac_labs <- sort(replicate(2,
                             paste(sample(letters, size = 7), collapse = "")))
y_fac_rlabs <- y_fac_labs[2:1]

# outcome variables:
ovars <- list(
  int = y1, # integer format
  rint = -y1 + 1, # reverse coded
  fac = factor(y1,
               levels = 0:1,
               labels = y_fac_labs), # factor (alphabetical)
  rfac = factor(y1,
                levels = 0:1,
                labels = y_fac_rlabs), # factor (non-alpha) - error expected.
  char = y_fac_labs[y1 + 1], # character (one way)
  rchar = y_fac_rlabs[y1 + 1] # charcter (the other)
)

# get all output:
res <- lapply(ovars, function(y) {
  try(parse_dCVnet_input(y = y, data = x1, family = "binomial"),
      silent = TRUE)
})

test_that("binomial y is formatted as expected", {
  expect_identical(res[[1]]$y, as.factor(ovars[[1]]))
  expect_identical(res[[2]]$y, as.factor(ovars[[2]]))
  expect_identical(res[[3]]$y, ovars[[3]])
  expect_identical(class(res[[4]]), "try-error")
  expect_identical(as.character(res[[5]]$y), ovars[[5]])
  expect_identical(as.character(res[[6]]$y), ovars[[6]])
})

# ~ multinomial -------------------------------------------------------------

# make some factor labels
set.seed(42)
m1 <- sample(1:5, size = 100, replace = TRUE)
m1.x <- matrix(0, nrow = length(m1), ncol = 5)
m_fac_labs <- sort(replicate(5,
                             paste(sample(letters, size = 7), collapse = "")))
m_fac_rlabs <- m_fac_labs[5:1]

# outcome variables:
ovars <- list(
  int = m1, # integer format
  rint = -m1 + 6, # reverse coded
  fac = factor(m1,
               levels = 1:5,
               labels = m_fac_labs), # factor (alphabetical)
  rfac = factor(m1,
                levels = 1:5,
                labels = m_fac_rlabs), # factor (non-alpha) - error expected.
  char = m_fac_labs[m1], # character (one wam)
  rchar = m_fac_rlabs[m1] # charcter (the other)
)

# get all output:
res <- lapply(ovars, function(y) {
  try(parse_dCVnet_input(y = y, data = m1.x, family = "binomial"),
      silent = TRUE)
})

test_that("binomial y is formatted as expected", {
  expect_identical(res[[1]]$y, as.factor(ovars[[1]]))
  expect_identical(res[[2]]$y, as.factor(ovars[[2]]))
  expect_identical(res[[3]]$y, ovars[[3]])
  expect_identical(class(res[[4]]), "try-error")
  expect_identical(as.character(res[[5]]$y), ovars[[5]])
  expect_identical(as.character(res[[6]]$y), ovars[[6]])
})
