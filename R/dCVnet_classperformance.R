

# classperformance S3 generic ---------------------------------------------


# a dCVnet S3 generic: classperformance
#   returns a merged dataframe or a list of data frames
#   where the dataframe contains raw performance information.


#  ~ making classperformance ----------------------------------------------

#' classperformance
#'
#' extracts a standardised classification performance table for a model.
#'
#' @name classperformance
#'
#' @param x an object from which binary class performance can be extracted.
#' @param ... arguments to pass on
#'
#' @return a classperformance object, is a dataframe (or list of dataframes)
#'     with the following columns:
#'    \itemize{
#'    \item{reference - the known 'true' class of the observation}
#'    \item{prediction - the model prediction for a case.
#'          for dCVnet this is the result of predict(model, type = "response")
#'          for "binary" probability of the non-positive class assigned
#'        by the model (used for AUROC)}
#'    \item{classification - for binomial and multinomial models this is
#'        the predicted class assigned by the model}
#'    \item{label - a grouping variable when predictions come from more than
#'        one source, e.g. multiple reps}
#'    }
#'
#' @export
classperformance <- function(x, ...) {
  UseMethod("classperformance", x)
}


#' classperformance.list
#' @describeIn classperformance classperformance for \code{\link{list}} object
#'     (treats it as a \code{\link{dCVnet}} object)
#' @export
classperformance.list <- function(x, ...) {
  # if we are given a list, lets assume it is a dCVnet style object
  #   such as dCVnet$final
  if (!"performance" %in% names(x)) {
    stop("not a suitable list for classperformance")
  }
  expected <- c("label", "reference", "prediction")
  if ( any(!expected %in% names(x$performance)) ) {
    stop("required columns for classperformance missing.")
  }
  classperformance.dCVnet(x, ...)
}


#' classperformance.classperformance
#' @describeIn classperformance returns self (allows list/dataframe conversion)
#' @export
classperformance.classperformance <- function(x, as.data.frame = TRUE, ...) {
  if ( as.data.frame && ! "data.frame" %in% class(x) ) {

    xfac <- as.factor(unlist(lapply(x, "[[", "label"), use.names = FALSE))
    return(structure(unsplit(x, xfac),
                     class = c("classperforamnce", "data.frame")))
  }
  if ( ! as.data.frame && "data.frame" %in% class(x) ) {
    x <- split(x, x$label)
    return(structure(x, class = c("classperformance", "list")))
  }
  return(x)
}


#' classperformance.dCVnet
#' @describeIn classperformance classperformance for \code{\link{dCVnet}} object
#' @param as.data.frame return a data.frame instead of a list of
#'     \code{\link{classperformance}} objects.
#' @export
classperformance.dCVnet <- function(x, as.data.frame = TRUE, ...) {
  if ( identical(as.data.frame, TRUE) ) {
    return(x$performance)
  } else {
    R <- split(x$performance, x$performance$label)
    return(structure(R, class = c("classperformance", "list")))
  }
}

#' classperformance.glm
#' @describeIn classperformance classperformance for \code{\link[stats]{glm}}
#'     object
#' @param label specify a label for the output
#' @param threshold for logistic regression use a threshold other than 0.5.
#' @export
classperformance.glm <- function(x,
                                 as.data.frame = TRUE,
                                 label = deparse(substitute(x)),
                                 threshold = 0.5, ...) {
  # Return (labelled) prediction dataframe from a glm
  #     given a threshold (default = 0.5):
  outcome <- as.character(x$terms[[2]])

  rwid <- rownames(x$data)

  lvl <- levels(x$data[[outcome]])
  classification <- as.numeric(stats::fitted(x) > threshold) + 1
  classification <- factor(lvl[classification], levels = lvl)

  reference <- x$data[[outcome]]
  prediction <- fitted(x)

  R <- data.frame(rowid = rwid,
                  reference = reference,
                  prediction = prediction,
                  classification = classification,
                  label = label)
  # return merged df or list.
  #   glms are always lists of length 1.
  if ( as.data.frame ) {
    return(structure(R, class = c("classperformance", "data.frame")))
  } else {
    return(structure(list(R), class = c("classperformance", "list")))
  }
}

#' classperformance.glm
#' @describeIn classperformance classperformance for glmlist from
#'     \code{\link{reflogreg}} object
#' @export
classperformance.glmlist <- function(x, as.data.frame = TRUE, ...) {
  # applies pobj.glm to a list of glms.

  class_list <- c("classperformance", "list")
  class_df <- c("classperformance", "data.frame")

  R <- lapply(seq_along(x), function(i) {
    # for a list we force return of a dataframe as we wrap in a list anyway.
    classperformance.glm(x[[i]],
                         as.data.frame = TRUE,
                         label = names(x)[i], ...)
  })
  names(R) <- names(x)
  if ( !as.data.frame ) return(structure(R, class = class_list))
  R <- as.data.frame(data.table::rbindlist(R))
  rownames(R) <- NULL
  return(structure(R, class = class_df))
}


#  ~ utilising class performance ----------------------------------------------

# Simple print function for class performance objects.
#' @export
print.classperformance <- function(x, ...) {
  if ( "list" %in% class(x) ) {
    type <- "list"
    n_models <- length(x)
    mod_labs <- names(x)
    px <- x[[1]]
  } else {
    type <- "data.frame"
    mod_labs <- unique(as.character(x$label))
    n_models <- length(mod_labs)
    sel <- as.character(x$label)
    sel <- (sel == sel[1])
    px <- subset(x, sel)
  }
  cat(paste("\nClassperformance object of type: ", type, "\n"))
  cat(paste("\tcontains results of", n_models, "model(s)\n"))
  cat(paste("\nOutcomes:"))
  print(table(px$reference))
  cat(paste("\nModels:\n\t"))
  cat(mod_labs)
  invisible(x)
}


#' summary.classperformance
#'
#' Calculates classification performance table and
#'     two-class classification metrics for a
#'     \code{\link{classperformance}} object.
#' @param object a \code{\link{classperformance}} object.
#' @param label a label can be assigned here.
#'      (Warning - setting a length 1 vector will concatenate multiple reps.)
#' @param ... additional arguments (not currently used)
#' @export
summary.classperformance <- function(object, label = NA, ...) {
  # Function assigns a label if asked to (for multi performance mode.)
  #   In multiperformance mode it uses label to produce columns of data.
  # If label is 'None' then this column is removed.

  if ( "list" %in% class(object) ) {
    # convert lists to data.frames.
    object <- structure(data.frame(do.call(rbind, object)),
                        class = c("classperformance", "data.frame"))
  }

  # Check structure:
  test <- names(object)
  test_cols <- c("reference", "classification", "label")
  test <- any(!(test_cols %in% test))
  if ( test ) {
    cat(names(object))
    stop("Check input")
  }

  if ( !is.na(label) ) object$label <- label

  # Two methods:
  .single_cpsummary <- function(performance) {
    # First categorical classification:
    A <- tidy_confusionmatrix(
      caret::confusionMatrix(
        data = performance$classification,
        reference = performance$reference))

    # Next add the AUC:
    B <- ModelMetrics::auc(actual = performance$reference,
                           predicted = performance$prediction)
    B <- pmax(B, 1 - B)
    B <- data.frame(Measure = "AUROC", Value = B)

    B$label <- A$label <- unique(performance$label)
    return(rbind(A, B))
  }

  .multi_cpsummary <- function(performance) {
    R <- lapply(seq_along(unique(performance$label)),
                function(i){
                  rr <- as.character(unique(performance$label)[i])
                  dd <- performance[performance$label == rr, ]
                  R <- summary.classperformance(dd, rr)
                  # Parse back to long.
                  R$label <- NULL # rr
                  names(R)[2] <- rr #"Value"
                  return(R)
                })
    R <- Reduce(function(x, y) merge(x, y, by = "Measure", sort = FALSE), R)
    return(R)
  }
  # If we lack a label col assign it and return
  #   (for multiclassperf functionality)
  # If there is a single label return single, otherwise return mulit.
  if ( is.null(object$label) ) {
    R <- .single_cpsummary(object)
  } else {
    if ( length(unique(object$label)) == 1 ) {
      R <- .single_cpsummary(object)
    } else {
      R <- .multi_cpsummary(object)
    }
  }
  # Option : remove label if it is there.
  if ( label %in% "None" ) R$label <- NULL
  return(R)
}

#' report_classperformance_summary
#'
#' extracts classperformance from a dCVnet object
#'     calculates classification statistics and
#'     provides a summary of
#'     \code{\link[base]{mean}},
#'     \code{\link[stats]{sd}},
#'     \code{\link[base:Extremes]{min}}, and \code{\link[base:Extremes]{max}}
#'
#' @name report_classperformance_summary
#'
#' @param dCVnet_object result from a call to \code{\link{dCVnet}}
#'
#' @return a data.frame of summarised and raw performance statistics.
#'
#' @export
report_classperformance_summary <- function(dCVnet_object) {

  outernreps <- length(unique(dCVnet_object$performance$label))

  if ( outernreps == 1 ) {
    # Simpler frame can be returned if only one outer loop:
    ols <- summary(classperformance(dCVnet_object), label = "None")
    names(ols)[2] <- "Rep1"
    return(ols)
  }

  ols <- summary(classperformance(dCVnet_object))

  summary_measures <- c("mean", "sd", "min", "max")
  names(summary_measures) <- summary_measures

  S <- lapply(summary_measures, function(M) {
    apply(ols[, -1, drop = FALSE], 1, function(i) {
      # suppress warning messages to quietly deal with NaNs in Mcnemar P-value.
      suppressWarnings(get(M)(i, na.rm = TRUE))
    })
  } )

  S <- do.call(data.frame, list(Measure = ols[, 1, drop = FALSE], S))

  S <- data.frame(S,
                  "..." = " - ",
                  ols[, -1, drop = FALSE],
                  stringsAsFactors = FALSE)

  return(S)
}



#' casesummary.classperformance
#'
#' What proportion of the time were subjects correctly classified in a
#'     \code{\link{classperformance}} object.
#' @param object a \code{\link{classperformance}} object?
#' @param type What should be returned?
#'                 \itemize{
#'                 \item{\code{data} - The true and estimated classifications.}
#'                 \item{\code{summary} - The mean proportion correct}
#'                 }
#' @param ... additional arguments (not currently used)
#' @export
casesummary.classperformance <- function(object,
                                         type = c("both",
                                                  "data",
                                                  "summary"),
                                         ...) {
  type <- match.arg(type)

  object <- as.data.frame(object)
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
  repdata <- lapply(labs, function(rep){
    object[object$label == rep, ]
  })
  names(repdata) <- names(labs)

  Rleft <- repdata[[1]][, c("rowid", "reference")]
  Rbits <- data.frame(lapply(repdata, function(kk) kk$classification),
                      stringsAsFactors = FALSE)
  R.data <- data.frame(lapply(Rbits, function(x) {
    as.numeric(x == Rleft$reference)
  } ))

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
