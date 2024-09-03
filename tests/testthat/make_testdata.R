
# s1: a single dCVnet model -----------------------------------------------

library(dCVnet)
make_dCVnet <- function() {
  data(BinomialExample, package = "glmnet")
  dCVnet(y = BinomialExample$y,    #nolint
         data = BinomialExample$x, #nolint
         family = "binomial",
         alphalist = 0.5,
         nrep_inner = 1)
}

m <- make_dCVnet()
save(m, file = "tests/testthat/data/model.RData", compress = TRUE)


# s2: examples from glmnet ------------------------------------------------

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
examples <- do.call(list, args = lapply(exlist, FUN = as.name))
examples <- setNames(examples, exlabs)

examples <- examples[!names(examples) %in% c("sparse")]

# convert cox y to Surv
examples$cox$y <- survival::Surv(time = examples$cox$y[, 1],
                                 event = examples$cox$y[, 2],
                                 type = "right")

save(examples,
     file = "tests/testthat/data/examples.RData",
     compress = TRUE)
