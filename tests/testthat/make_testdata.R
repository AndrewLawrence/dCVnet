
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
