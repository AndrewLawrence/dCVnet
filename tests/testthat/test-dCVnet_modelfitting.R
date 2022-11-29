
# Setup:

# glmnet 4.1.4 changed how examples work, now the example creates a
#   variable of the name of the example which is a two element list
#   containing x and y.

# identify available glmnet examples:
exlist <- data(package = "glmnet")$results[, "Item"]
exlist <- exlist[grepl("Example", exlist)]

# load the data for each:
for (i in exlist) {
  data(package = "glmnet", list = i)
}

exlabs <- tolower(gsub("Example", "", exlist))
exlabs[exlabs == "quickstart"] <- "gaussian"
exlabs[exlabs == "multigaussian"] <- "mgaussian"


# merge
examples <- do.call(list, args = sapply(exlist, as.name))
examples <- setNames(examples, exlabs)

examples <- examples[!names(examples) %in% c("sparse")]

# convert cox y to Surv
examples$cox$y <- survival::Surv(time = examples$cox$y[, 1],
                                 event = examples$cox$y[, 2],
                                 type = "right")


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
test_that(
  "multialpha objects pass test", {
  for ( i in seq_along(ma.glmnets) ) {
    expect_true(is_multialpha.repeated.cv.glmnet(ma.glmnets[[i]]))
  }
  }
)
