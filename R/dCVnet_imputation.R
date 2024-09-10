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
  .pp_fit_mean <- function(x, family, y = NULL) {
    if (!is.null(y)) {
      x <- impy_dat_merger(x = x, y = y, family = family)
    }
    caret::preProcess(x, method = c("center", "scale"))
  }

  .pp_apply_mean <- function(x, newdata, family, newy = NULL) {
    if ( is.null(newy) ) {
      r <- as.matrix(predict(x, newdata = newdata))
      r[is.na(r)] <- 0.0
      return(r)
    }
    new_comb <- impy_dat_merger(x = newdata, y = newy, family = family)
    r <- predict(x, new_comb)
    r[is.na(r)] <- 0.0

    attr(r, "family") <- attr(new_comb, "family")
    if ( family == "cox" ) attr(r, "Survtype") <- attr(new_comb, "Survtype")
    attr(r, "ny") <- attr(new_comb, "ny")

    # undo centering / scaling for y:
    for ( i in seq.int(attr(r, "ny")) ) {
      if ( ! is.factor(r[, i]) ) {
        r[, i] <- (r[, i] * x$std[[i]]) + x$mean[[i]]
      }
    }

    r <- impy_dat_unmerger(r)
    r$x
  }



  # knn imputation:
  .pp_fit_caretknn <- function(x, family, y = NULL) {
    if (!is.null(y)) {
      x <- impy_dat_merger(x = x, y = y, family = family)
    }
    caret::preProcess(x, method = c("center", "scale", "knnImpute"))
  }

  .pp_apply_caret <- function(x, newdata, family, newy = NULL) {
    if ( is.null(newy)) {
      return(as.matrix(predict(x, newdata = newdata)))
    }
    new_comb <- impy_dat_merger(x = newdata, y = newy, family = family)
    r <- predict(x, new_comb)

    attr(r, "family") <- attr(new_comb, "family")
    if ( family == "cox" ) attr(r, "Survtype") <- attr(new_comb, "Survtype")
    attr(r, "ny") <- attr(new_comb, "ny")

    # undo centering / scaling for y:
    for ( i in seq.int(attr(r, "ny")) ) {
      if ( ! is.factor(r[, i]) ) {
        r[, i] <- (r[, i] * x$std[[i]]) + x$mean[[i]]
      }
    }

    r <- impy_dat_unmerger(r)
    r$x
  }

  # missForestPredict imputation:
  .pp_fit_mfp <- function(x, family, y = NULL) {
    requireNamespace("missForestPredict", quietly = TRUE)
    if ( ! is.null(y) ) {
      x <- impy_dat_merger(x = x, y = y, family = family)
    }

    mfp <- missForestPredict::missForest(as.data.frame(x),
                                         save_models = TRUE, verbose = FALSE)
    PPx <- caret::preProcess(
      missForestPredict::missForestPredict(mfp,
                                           newdata = as.data.frame(x)),
      method = c("center", "scale")
    )
    return(list(missForest = mfp, PPx = PPx))
  }

  .pp_apply_mfp <- function(x, newdata, family, newy = NULL) {

    if ( is.null(newy) ) {

      newdata <- missForestPredict::missForestPredict(
        x[["missForest"]],
        newdata = as.data.frame(newdata)
      )

      return(as.matrix(predict(x[["PPx"]], newdata = newdata)))
    }

    new_comb <- impy_dat_merger(x = newdata,
                                y = newy,
                                family = family)

    newdata <- missForestPredict::missForestPredict(
      x[["missForest"]],
      newdata = as.data.frame(new_comb)
    )
    r <- predict(x[["PPx"]], newdata = newdata)

    attr(r, "family") <- attr(new_comb, "family")
    if ( family == "cox" ) attr(r, "Survtype") <- attr(new_comb, "Survtype")
    attr(r, "ny") <- attr(new_comb, "ny")

    # undo centering / scaling for y:
    for ( i in seq.int(attr(r, "ny")) ) {
      if ( ! is.factor(r[, i]) ) {
        r[, i] <- (r[, i] * x[["PPx"]]$std[[i]]) + x[["PPx"]]$mean[[i]]
      }
    }

    r <- impy_dat_unmerger(r)
    r[["x"]]
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


# In order to use y-variables in the imputation we need to
#   merge and unmerge y- and x- input without mangling
#   or loss of information.
# The following two functions should achieve this.
impy_dat_merger <- function(x, y, family) {
  ny <- NCOL(y)
  if ( family == "cox" ) {
    r <- data.frame(y = as.matrix(y), x)
    attr(r, "Survtype") <- attr(y, "type")
  } else {
    r <- data.frame(y = y, x)
    attr(r, "Survtype") <- NULL
  }
  attr(r, "family") <- family
  attr(r, "ny") <- ny
  r
}

impy_dat_unmerger <- function(x) {
  family <- attr(x, "family")
  Survtype <- attr(x, "Survtype")
  ny <- attr(x, "ny")

  if ( ny > 1 ) {
    ysel <- seq.int(NCOL(x)) %in% seq.int(ny)
    y <- as.matrix(x[, ysel])
    rownames(y) <- NULL
    colnames(y) <- gsub("^y.", "", colnames(y))

    x <- as.matrix(x[, !ysel])
  } else {
    y <- x[, 1]
    x <- as.matrix(x[, -1])
  }

  if ( family == "cox" ) {
    if ( ny == 2 ) {
      y <- survival::Surv(time = as.vector(y[, 1]),
                          event = as.vector(y[, 2]),
                          type = Survtype)
    } else {
      y <- survival::Surv(time = as.vector(y[, 1]),
                          time2 = as.vector(y[, 2]),
                          event = as.vector(y[, 3]),
                          type = Survtype)
    }

  }
  list(x = x,
       y = y)
}
