# Functions related to the implementation of imputation in dCVnet

# The preproc_imp_functions function returns functions.
#   It produces a list of two functions selected appropriate to
#     the opt.imputation_method argument.
#   The functions are:
#       fit - for running preprocessing + imputation
#             the "fit" function returns an object.
#       apply - for making predictions using the object
#               returned from fit and some newdata.
# preproc_imp_functions does not use y in imputation or prediction.
preproc_imp_functions <- function(opt.imputation_method) {
  # Mean imputation:
  .pp_fit_mean <- function(x) {
    caret::preProcess(x, method = c("center", "scale"))
  }
  .pp_apply_mean <- function(x, newdata) {
    newdata[is.na(newdata)] <- 0.0
    as.matrix(predict(x, newdata = newdata))
  }
  # knn imputation:
  .pp_fit_caretknn <- function(x) {
    caret::preProcess(x, method = c("center", "scale", "knnImpute"))
  }
  .pp_apply_caret <- function(x, newdata) {
    as.matrix(predict(x, newdata = newdata))
  }
  # missForestPredict imputation:
  .pp_fit_mfp <- function(x) {
    requireNamespace("missForestPredict", quietly = TRUE)
    mfp <- missForestPredict::missForest(as.data.frame(x),
                                         save_models = TRUE, verbose = FALSE)
    PPx <- caret::preProcess(
      missForestPredict::missForestPredict(mfp,
                                           newdata = as.data.frame(x)),
      method = c("center", "scale")
    )
    list(missForest = mfp, PPx = PPx)
  }
  .pp_apply_mfp <- function(x, newdata) {
    newdata <- missForestPredict::missForestPredict(
      x[["missForest"]],
      newdata = as.data.frame(newdata)
    )
    as.matrix(predict(x[["PPx"]], newdata = newdata))
  }
  # Main:
  pp_fit <- switch(
    opt.imputation_method,
    mean = .pp_fit_mean,
    knn = .pp_fit_caretknn,
    missForestPredict = .pp_fit_mfp
  )
  pp_apply <- switch(
    opt.imputation_method,
    mean = .pp_apply_mean,
    knn = .pp_apply_caret,
    missForestPredict = .pp_apply_mfp
  )
  return(list(fit = pp_fit, apply = pp_apply))
}

