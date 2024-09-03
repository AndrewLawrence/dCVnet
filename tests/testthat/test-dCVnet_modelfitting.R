
load(file = "data/examples.RData")

# Run dCVnet parsing for examples:
parsed <- mapply(function(x, y) {
  parse_dCVnet_input(data = x$x,
                     y = x$y,
                     family = y)
},
examples,
names(examples),
SIMPLIFY = FALSE)

# Run multi-alpha.repeated.cv.glmnets for the examples:
capture.output(
  ma.glmnets <- lapply(parsed, function(xxx) {
    suppressWarnings(do.call(multialpha.repeated.cv.glmnet,
                             args = list(x = xxx$x_mat,
                                         y = xxx$y,
                                         family = xxx$family,
                                         opt.keep_models = "best",
                                         k = 3,
                                         nrep = 1,
                                         alphalist = c(0.2, 0.4, 0.6))))
  } )
)

# Tests:
test_that("multialpha objects pass test", {
  for (i in seq_along(ma.glmnets)) {
    expect_true(is_multialpha.repeated.cv.glmnet(ma.glmnets[[i]]))
  }
})
