

# performance S3 generic ---------------------------------------------

# a dCVnet S3 generic: performance
#   returns a merged dataframe or a list of data frames
#   where the dataframe contains raw performance information:
#     - reference       (actual 'true' outcome)
#     - prediction      (predicted outcome,
#                         will be probability for classifications)
#     - classification  (the classification derived from prediction,
#                         if applicable)
#     - label           (a grouping label e.g. rep10fold3, model1)
#     - rowid           (observation labels from the rowids of the data)
#
#   The attributes of the performance object must contain the model
#     family
#
# Conceptually a performance object contains all the data required to
#   evaluate the performance of a model, or models.
#
#   The actual performance measures are then calculated by generic
#       print / summary methods.
#
#   Thus, to construct a performance object you need all the input to the
#     predict function appropriate for that model plus the actual outcome
#     being predicted.
#
#   Further (optional) extensions are:
#     labels  - allowing performance from different models to be held
#               in a single object.
#     rowid   - allowing performance to be broken down by subject, and
#               simplifying the process of merging / augmentation of raw data.
#
# Practically the main difference between a performance object and a
#   predict object is the performance object is tidy, and must also encode
#   the observed outcomes associated with newx - predict methods only require
#   newx. We also add rownames / caselabels when these are present.
#
# important notes:
#   these functions are in development for model families other than binomial
#   y formats vary substantially between the model families
#   notably:
#     cox has both outcome (status) and censoring time in y, but
#       results of predict are a vector of relative risks
#     multinomial can take either a k>2 factor or a matrix with k columns
#       and returns a matrix with k-columns
#     mgaussian takes a k-column matrix and returns a k-column matrix.
#
# tidy_predict.glmnet does most of the heavy lifting here for
#   a glmnet objects
#   performance should be recoded to make more use of this.

#  ~ making performance ----------------------------------------------

#' performance
#'
#' extracts a standardised classification performance table for a model
#'     or list of models.
#'
#' @name performance
#'
#' @param x an object from which prediction performance can be extracted.
#' @param ... arguments to pass on
#'
#' @return a performance object, is a dataframe (or list of dataframes)
#'     with the following columns:
#'    \itemize{
#'    \item{reference - the known 'true' class of the observation}
#'    \item{prediction - the model prediction for a case.
#'          for dCVnet this is the result of predict(model, type = "response")
#'          With "binary" response predictions are the predicted probability
#'          of the non-reference class (used for AUROC)}
#'    \item{classification - for binomial and multinomial models this is
#'        the predicted class assigned by the model}
#'    \item{label - a grouping variable when predictions come from more than
#'        one source, e.g. multiple reps}
#'    }
#'
#' @export
performance <- function(x, ...) {
  UseMethod("performance")
}

#' performance.default
#'
#' @rdname performance
#' @description Default function behaviour assumes input is a
#'     \code{\link{dCVnet}} object (and fails if this is not reasonable).
#' @export
performance.default <- function(x, ...) {
  # if we are given a list, lets assume it is a dCVnet style object
  #   such as dCVnet$final

  # Because we're not certain then check the structure:
  if (!"performance" %in% names(x)) {
    stop("not a suitable list for performance (wrong names)")
  }
  if ( inherits(try( family(x), "try-error")) ) {
    stop("model family missing")
  }
  expected <- c("label", "reference", "prediction")
  if ( any(!expected %in% names(x$performance)) ) {
    stop("required columns for performance missing.")
  }

  performance.dCVnet(x, ...)
}



#' performance.dCVnet
#'
#' @rdname performance
#' @param as.data.frame return a data.frame instead of a list of
#'     \code{\link{performance}} objects.
#' @export
performance.dCVnet <- function(x, as.data.frame = TRUE, ...) {
  if ( identical(as.data.frame, TRUE) ) {
    R <- x$performance
  } else {
    R <- structure(split(x$performance, x$performance$label),
                   class = c("performance", "list"))
  }
  attr(R, which = "family") <- family(x)
  return(R)
}


#' performance.performance
#'
#' @rdname performance
#' @description Applying performance to a performance object
#'     allows conversion between list/dataframe format.
#' @export
performance.performance <- function(x, as.data.frame = TRUE, ...) {
  R <- x # fall-back
  if ( as.data.frame && !inherits(x, "data.frame") ) {

    xfac <- as.factor(unlist(lapply(x, "[[", "label"), use.names = FALSE))
    R <- structure(unsplit(x, xfac),
                   class = c("performance", "data.frame"))
  }
  if ( ! as.data.frame && inherits(x, "data.frame") ) {
    R <- structure(split(x, x$label),
                   class = c("performance", "list"))
  }
  attr(R, which = "family") <- family(x)
  return(R)
}

#' performance.glm
#'
#' For glm objects wraps \link{predict.glm} if newdata is specified.
#'
#' @rdname performance
#' @param label specify a label for the output
#' @param newdata evaluate performance in new data
#' @param threshold for binomial family glm - use specified threshold
#'     for predicting class from probability.
#' @export
performance.glm <- function(x,
                            as.data.frame = TRUE,
                            label = deparse(substitute(x)),
                            threshold = 0.5,
                            newdata = NULL,
                            ...) {
  # dataframe of prediction results from a glm
  #     given a threshold (default = 0.5).

  # extract model.frame:
  mf <- stats::model.frame(x)
  mf.y <- stats::model.response(mf)
  # name of outcome (always leftmost in model.frame):
  outcome <- colnames(mf)[[1]]
  # name of family
  familyname <- family(x)$family

  lvl <- NULL
  if ( is.factor(mf.y) ) {
    lvl <- levels(mf.y) # NULL if not factor
  }

  if ( is.null(newdata) ) {
    # extraction from model object:
    rwid <- rownames(mf)
    prediction <- stats::fitted(x)
    reference <- mf.y
  } else {
    ## There will be problems if y in the data doesn't have the same class
    ##        and levels as y in the newdata:
    stopifnot(identical(class(mf.y),
                        class(newdata[[outcome]])))
    if ( ! is.null(lvl) ) {
      stopifnot(identical(levels(mf.y),
                          levels(newdata[[outcome]])))
    }

    rwid <- rownames(newdata)
    if ( is.null(rwid) ) rwid <- seq_len(NROW(newdata))
    prediction <- predict(x,
                          newdata = newdata,
                          type = "response", ...)
    reference <- newdata[[outcome]]
  }
  classification <- rep(NA, length(mf.y))

  if ( familyname %in% c("binomial") && length(unique(mf.y)) < 3 ) {
    # convert probabilities to classifications:
    classification <- as.integer(prediction > threshold)
    if ( !is.null(lvl) ) {
      # use factor levels if they exist:
      classification <- factor(lvl[classification + 1L], levels = lvl)
    } else {
      # if we don't have levels we have glm data which *must* be 0 1:
      classification <- factor(classification, levels = c(0, 1))
      # also cooerce the reference data:
      reference <- factor(reference, levels = c(0, 1))
    }
  }

  R <- structure(data.frame(rowid = rwid,
                            reference = reference,
                            prediction = prediction,
                            classification = classification,
                            label = label,
                            stringsAsFactors = FALSE),
                 class = c("performance", "data.frame"),
                 family = familyname)
  # return merged df or list.
  return(performance(R, as.data.frame = as.data.frame))
}


#' performance.glmlist
#'
#' @rdname performance
#' @export
performance.glmlist <- function(x, as.data.frame = TRUE, ...) {
  # applies pobj.glm to a list of glms.

  class_list <- c("performance", "list")
  class_df <- c("performance", "data.frame")

  familyname <- family(x[[1]])$family

  R <- lapply(seq_along(x), function(i) {
    performance.glm(x[[i]],
                    as.data.frame = TRUE,
                    label = names(x)[i], ...)
  })
  names(R) <- names(x)
  if ( !as.data.frame ) return(structure(R, class = class_list))
  R <- as.data.frame(data.table::rbindlist(R), stringsAsFactors = FALSE)
  rownames(R) <- NULL
  return(structure(R, class = class_df, family = familyname))
}


# WIP: write a performance.glm_net_ that calls to tidy_predict.glmnet.
#
#

#  ~ utilising class performance ----------------------------------------------

# Simple print function for class performance objects.
#' @export
print.performance <- function(x, ...) {
  familyname <- family(x)
  if ( inherits(x, "list") ) {
    type <- "list"
    n_models <- length(x)
    mod_labs <- names(x)
  } else if (inherits(x, "data.frame")) {
    type <- "data.frame"
    mod_labs <- unique(as.character(x$label))
    n_models <- length(mod_labs)
  } else {
    stop("object must inherit list or data.frame")
  }
  cat(paste("\nperformance object of type: ", type, "\n"))
  cat(paste("\tcontains results of",
            n_models, "reps of", familyname, "family model(s)\n"))
  cat(paste("\nOutcomes:\n"))
  print(describe_y_from_performance(x))
  cat(paste("\nModels:\n\t"))
  cat(mod_labs)
  invisible(x)
}

#' @name InternalPerformanceSummaryFunctions
#'
#' @title Internal summary functions for performance objects from
#'     different model families.
#'
#' @rdname InternalPerformanceSummaryFunctions
NULL

#' @describeIn InternalPerformanceSummaryFunctions
#'     Used for binomial and multinomial families
perf_nomial <- function(object,
                        short = FALSE,
                        pvprevalence = "observed") {
  if (identical(pvprevalence, "observed")) {
    pvprevalence <- NULL
  }

  # First categorical classification:
  A <- tidy_confusionmatrix(
    caret::confusionMatrix(
      data = object$classification,
      reference = object$reference,
      # ~~~~~~~~~~~~~
      # Factor levels
      # ~~~~~~~~~~~~~
      # The ordering of factor levels is important for calculating
      #   sensitivity, specificity, PPV, NPV etc.
      # caret by default (and confusion matrices in general) follows the
      #   convention that the class being predicted (e.g. diseased subjects)
      #   is stored in the first level of the factor.
      #   However glm / model.matrix uses treatment-coding: the
      #   first level of y is the reference class (e.g. control subjects).
      #   This behaviour is helpful when interpreting model coefficients
      #   as they represent the deviation from the reference such that
      #   a positive coefficient indicates an increased probability of
      #   the 'positive' class.
      # glmnet: to add to the confusion glmnet will reorder binary factors
      #   such that levels are alphabetical.
      # dCVnet coerces input data into an alphabetical factor with the
      #   reference class as the first level, so when calculating
      #   classification performance we must force level 2 as the
      #   'positive' class (where positive is taken in the sense of
      #   'positive predictive value'). This is done here:
      positive = levels(object$reference)[[2]],
      prevalence = pvprevalence
    )
  )

  # Next add the AUC:
  if (length(unique(object$reference)) > 2) {
    mcols <- colnames(object)[which(grepl("^prediction",
                                          colnames(object)))]
    mlvl <- gsub("prediction", "", mcols)
    B <- ModelMetrics::mauc(actual = object$reference,
                            predicted = object[, mcols])
    B <- data.frame(Measure = c("multiclass OvR AUROC",
                                paste("OvR AUROC", mlvl)),
                    Value = unlist(B))
  } else {
    B <- ModelMetrics::auc(actual = object$reference,
                           predicted = object$prediction)
    B <- data.frame(Measure = "AUROC",
                    Value = B)
  }
  # and the Brier Score:
  if (length(unique(object$reference)) < 3) {
    if (is.numeric(object$reference)) {
      Bs <- ModelMetrics::brier(actual = object$reference,
                                predicted = object$prediction)
    } else {
      # for factors need conversion:
      Bs <-
        ModelMetrics::brier(actual = as.numeric(object$reference) - 1,
                            predicted = object$prediction)
    }
    Bs <- data.frame(Measure = "Brier", Value = Bs)
    B <- rbind(B, Bs)
  }

  B$label <- A$label <- unique(object$label)
  R <- rbind(A, B)

  if ( family(object) == "binomial" ) {
    min_vars <- c("Accuracy", "Sensitivity",
                  "Specificity", "Balanced Accuracy",
                  "AUROC")
  } else {
    lvls <- levels(object$reference)
    min_vars <- c("Accuracy",
                  paste("Class:", lvls, "Balanced Accuracy"),
                  "multiclass OvR AUROC",
                  paste("OvR AUROC", lvls))
  }

  if ( short ) {
    R <- R[R$Measure %in% min_vars, ]
  }
  return(R)
}

#' @describeIn InternalPerformanceSummaryFunctions
#'     Used for gaussian and poisson families
#' @importFrom ModelMetrics rmse mae rmsle
perf_cont <- function(object,
                      short = FALSE,
                      pvprevalence = "observed") {

  f <- family(object)

  ysd <- sd(object$reference)

  cval <- cor(object$reference,
              object$prediction)

  mod <- coef(lm(reference ~ prediction, data = as.data.frame(object)))

  # Using ModelMetrics:
  R <- c(
    RMSE = ModelMetrics::rmse(
      actual = object$reference,
      predicted = object$prediction
    ),
    MAE = ModelMetrics::mae(
      actual = object$reference,
      predicted = object$prediction
    ),
    r = cval,
    r2 = cval ^ 2,
    cal_Intercept = mod[[1]],
    cal_Slope = mod[[2]]
  )

  # Add brier score:
  R["Brier"] <- R["RMSE"]^2

  # Add scaled versions of RMSE & MAE:
  R["SDScaledRMSE"] <- R["RMSE"] / ysd
  R["SDScaledMAE"] <- R["MAE"] / ysd

  if ( f == "poisson" ) {
    log_dat <- data.frame(reference = log1p(object$reference),
                          prediction = log1p(object$prediction))

    log_cval <- cor(log_dat$reference, log_dat$prediction)
    log_mod <- coef(lm(reference ~ prediction, data = log_dat))

    R <- append(R, c(logged_r = log_cval,
                     logged_r2 = log_cval^2,
                     RMSLE = ModelMetrics::rmsle(actual = object$reference,
                                                 predicted = object$prediction),
                     logged_cal_Intercept = log_mod[[1]],
                     logged_cal_Slope = log_mod[[2]]))
  }

  R <- data.frame(Measure = names(R), Value = R, label = object$label[[1]])
  row.names(R) <- NULL

  R
}

#' @describeIn InternalPerformanceSummaryFunctions
#'     Used for cox models
#' @importFrom Hmisc rcorr.cens
#' @importFrom survival Surv is.Surv
perf_cox <- function(object,
                     short = FALSE,
                     pvprevalence = "observed") {
  R <- Hmisc::rcorr.cens(
    x = object$prediction,
    S = survival::Surv(object$reference.Time,
                       object$reference.Status)
  )[1:3]

  R <-
    data.frame(Measure = names(R),
               Value = R,
               label = object$label[[1]])
  row.names(R) <- NULL

  R
}




#' @describeIn InternalPerformanceSummaryFunctions
#'     Used for multivariate gaussian models
#'
perf_mgaussian <- function(object,
                           short = FALSE,
                           pvprevalence = "observed") {
  # predictions to iterate through:
  cn <- colnames(object)

  pcols <- cn[which(grepl("prediction", cn))]
  rcols <- gsub("prediction", "reference", pcols)

  vlabs <- gsub("prediction", "", pcols)

  keep_cols <- cn[!cn %in% c(rcols, pcols)]

  R <- lapply(seq_len(length(pcols)),
              function(i) {
                obj <- object[, c(pcols[i], rcols[i], keep_cols)]
                colnames(obj) <-
                  c("prediction", "reference", keep_cols)
                attr(obj, "family") <- "gaussian"
                VAL <- perf_cont(obj)
                VAL$label <- pcols[i]
                VAL
              })

  RR <- R[[1]]
  RR$label <- NULL
  RR$Value <-
    (Reduce("+", x = lapply(R, function(x)
      x$Value))) / length(pcols)

  RR$Measure <- paste0("mean ", RR$Measure)

  R <- lapply(seq_len(length(pcols)),
              function(i) {
                x <- R[[i]]
                x$Measure <- paste0(vlabs[[i]], " ", x$Measure)
                x$label <- NULL
                return(x)
              })

  R <- list(RR,
            do.call("rbind", R))

  R <- do.call("rbind", R)

  if (short) {
    R <- R[grepl("RMSE", R$Measure), ]
  }

  return(R)
}


#' summary.performance
#'
#' Calculates classification performance table and
#'     two-class classification metrics for a
#'     \code{\link{performance}} object.
#'
#' @param object a \code{\link{performance}} object.
#' @param label a label can be assigned here.
#'      (Warning - setting a length 1 vector will concatenate multiple reps.)
#' @param short (bool) return a core set of performance measures.
#' @param pvprevalence argument for adjustment of PPV/NPV calculation for
#'     binomial or multinomial families. Either "observed" to use the observed
#'     prevalence, or a number \code{[0,1]} (for binomial),
#'     or a vector of length n_categories (for multinomial) containing numerics
#'     in \code{[0,1]}.
#' @param ... additional arguments (ignored)
#'
#' @export
summary.performance <- function(object,
                                label = NA,
                                short = FALSE,
                                pvprevalence = "observed",
                                ...) {
  # Function assigns a label if asked to (for multi performance mode.)
  #   In multiperformance mode it uses label to produce columns of data.
  # If label is 'None' then this column is removed.

  if ( "list" %in% class(object) ) {
    # convert lists to data.frames.
    object <- structure(data.frame(do.call(rbind, object),
                                   stringsAsFactors = FALSE),
                        class = c("performance", "data.frame"))
  }

  if ( !is.na(label) ) object$label <- label

  # Two methods:
  .single_cpsummary <- function(performance,
                                short = short,
                                pvprevalence = "observed") {
    short
    f <- family(performance)
    fxn <- switch(f,
                  binomial = perf_nomial,
                  multinomial = perf_nomial,
                  gaussian = perf_cont,
                  poisson = perf_cont,
                  cox = perf_cox,
                  mgaussian = perf_mgaussian,
                  stop("family not supported"))
    R <- fxn(performance,
               short = short,
               pvprevalence = pvprevalence)
    rownames(R) <- NULL
    return(R)
  }

  .multi_cpsummary <- function(performance,
                               short,
                               pvprevalence) {
    R <- lapply(seq_along(unique(performance$label)),
                function(i) {
                  rr <- as.character(unique(performance$label)[i])
                  dd <- performance[performance$label == rr, ]
                  R <- summary.performance(object = dd,
                                           label = rr,
                                           short = short,
                                           pvprevalence = pvprevalence)
                  # Parse back to long.
                  R$label <- NULL # rr
                  names(R)[2] <- rr #"Value"
                  return(R)
                })
    R <- Reduce(function(x, y) merge(x, y, by = "Measure", sort = FALSE), R)
    return(R)
  }

  # Processing:
  #   If we lack a label col assign it and return
  #     (for multiclassperf functionality)
  #   If there is a single label return single, otherwise return mulit.
  if ( is.null(object$label) ) {
    R <- .single_cpsummary(object, short = short, pvprevalence = pvprevalence)
  } else {
    if ( length(unique(object$label)) == 1 ) {
      R <- .single_cpsummary(object, short = short, pvprevalence = pvprevalence)
    } else {
      R <- .multi_cpsummary(object, short = short, pvprevalence = pvprevalence)
    }
  }
  # Option : remove label if it is there.
  if ( label %in% "None" ) R$label <- NULL
  return(R)
}


#' @describeIn family.dCVnet family for model performance objects
#' @export
family.performance <- function(object, ...) {
  attr(object, which = "family")
}


#' report_performance_summary
#'
#' extracts performance from a dCVnet object
#'     calculates classification statistics and
#'     provides a summary of
#'     \code{\link[base]{mean}},
#'     \code{\link[stats]{sd}},
#'     \code{\link[base:Extremes]{min}}, and \code{\link[base:Extremes]{max}}
#'
#' @name report_performance_summary
#'
#' @param dCVnet_object result from a call to \code{\link{dCVnet}}
#' @param short (bool) return a core set of performance measures.
#' @param pvprevalence allows calculation of PPV/NPV at different prevalences.
#'      set to "observed" to use the prevalence of the data.
#'      For binomial data use a single value, for multinomial use a
#'      named vector of prevalences with names as per the levels of y.
#'      Note: does not affect the presented prevalence value in the table.
#'
#' @return a data.frame of summarised and raw performance statistics.
#'
#' @importFrom psych fisherz fisherz2r
#'
#' @export
report_performance_summary <- function(dCVnet_object,
                                       short = FALSE,
                                       pvprevalence = "observed") {

  if ( inherits(dCVnet_object, "dCVnet") ) {
    outernreps <- length(unique(dCVnet_object$performance$label))
  } else {
    # assume it's a performance object and try to get label:
    outernreps <- length(unique(dCVnet_object$label))
  }

  if ( outernreps == 1 ) {
    # Simpler frame can be returned if only one outer loop:
    ols <- summary(performance(dCVnet_object),
                   label = "None",
                   short = short,
                   pvprevalence = pvprevalence)
    names(ols)[2] <- "Rep1"
    return(ols)
  }

  # extract performance measures for each label:
  ols <- summary(performance(dCVnet_object),
                 short = short,
                 pvprevalence = pvprevalence)

  fisher_measures <- c("r", "r2", "logged_r", "logged_r2")
  doing_fisherr2z <- any(fisher_measures %in% as.character(ols$Measure))

  # Fisher R to z business:
  if ( doing_fisherr2z ) {
    fish_sel <- as.character(ols$Measure) %in% fisher_measures
    ols[fish_sel, -1] <- psych::fisherz(ols[fish_sel, -1])
  }

  summary_measures <- c("mean", "sd", "min", "max")
  names(summary_measures) <- summary_measures

  S <- lapply(summary_measures, function(M) {
    apply(ols[, -1, drop = FALSE], 1, function(i) {
      # suppress warning messages to quietly deal with NaNs in Mcnemar P-value.
      suppressWarnings(get(M)(i, na.rm = TRUE))
    })
  } )


  S <- do.call(data.frame, list(Measure = ols[, 1, drop = FALSE],
                                S,
                                stringsAsFactors = FALSE))

  # Fisher R to z back transform:
  if ( doing_fisherr2z ) {
    # back transform the raw results:
    ols[fish_sel, -1] <- psych::fisherz2r(ols[fish_sel, -1])

    # and the averaged results:
    fish_sel <- as.character(S$Measure) %in% fisher_measures
    S[fish_sel, -1] <- psych::fisherz2r(S[fish_sel, -1])
  }

  S <- data.frame(S,
                  "..." = " - ",
                  ols[, -1, drop = FALSE],
                  stringsAsFactors = FALSE)

  if ( short ) {
    S <- S[, c("Measure", "mean", "sd", "min", "max")]
  }

  return(S)
}


#' casesummary.performance
#'
#' What proportion of the time were subjects correctly classified in a
#'     \code{\link{performance}} object.
#' @param object a \code{\link{performance}} object?
#' @param type What should be returned?
#'                 \itemize{
#'                 \item{\code{data} - The true and estimated classifications.}
#'                 \item{\code{summary} - The mean proportion correct}
#'                 }
#' @param ... additional arguments (not currently used)
#' @export
casesummary.performance <- function(object,
                                    type = c("both",
                                             "data",
                                             "summary"),
                                    ...) {
  type <- match.arg(type)

  object <- as.data.frame(object, stringsAsFactors = FALSE)
  labs <- unique(object$label)
  names(labs) <- make.names(labs, unique = TRUE)

  # sort the object by rowid (with numerical order if appropriate):
  if ( all(!is.na(suppressWarnings(as.numeric(object$rowid)))) ) {
    # if we can convert to numeric without any NAs, then force numeric:
    object <- object[order(as.numeric(object$rowid)), ]
  } else {
    # sort by character:
    object <- object[order(object$rowid), ]
  }

  # iterate over the labels (reps) to make a wide dataframe.
  repdata <- lapply(labs, function(rep) {
    object[object$label == rep, ]
  })
  names(repdata) <- names(labs)

  Rleft <- repdata[[1]][, c("rowid", "reference")]
  Rbits <- data.frame(lapply(repdata, function(kk) kk$classification),
                      stringsAsFactors = FALSE)
  R.data <- data.frame(lapply(Rbits, function(x) {
    as.numeric(x == Rleft$reference)
  } ), stringsAsFactors = FALSE)

  R <- switch(type,
              data = list(Rleft,
                          Rbits,
                          stringsAsFactors = TRUE),
              summary = list(Rleft,
                             prop_correct = rowMeans(R.data),
                             stringsAsFactors = TRUE),
              both = list(Rleft,
                          prop_correct = rowMeans(R.data),
                          `...` = "-",
                          Rbits,
                          stringsAsFactors = TRUE)
  )

  return(do.call(data.frame, R))
}

#' get_y_from_performance
#'
#' Extracts reference (y) values from a performance object
#'
#' @param object a performance object
#' @keywords internal
#' @noRd
get_y_from_performance <- function(object) {
  if ( family(object) == "cox"  ) {
    return(survival::Surv(
      time = object$reference.Time,
      event = object$reference.Status
    ))
  }
  # convert to "list" type and take first element:
  indata <- performance(object, as.data.frame = FALSE)[[1]]
  sel <- grepl("^reference", colnames(indata))
  y <- indata[, sel, drop = FALSE]
  class(y) <- "data.frame"
  return(y)
}


#' describe_y_from_performance
#'
#' Brief descriptives of y
#'
#' @param object a performance object
#' @keywords internal
#' @noRd
describe_y_from_performance <- function(object) {
  familyname <- family(object)
  y <- get_y_from_performance(object)
  describe_outcome(y, family = familyname)
}
