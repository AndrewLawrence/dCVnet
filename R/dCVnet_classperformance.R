

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
#'    \item{classification - the predicted class assigned by the model}
#'    \item{probability - the probability of the non-positive class assigned
#'        by the model (used for AUROC)}
#'    \item{label - a grouping variable when predictions come from more than
#'        one source, e.g. multiple reps}
#'    }
#'
#' @export
classperformance <- function(x, ...) {
  UseMethod("classperformance", x)
}

#' @export
classperformance.list <- function(x, ...) {
  # if we are given a list, lets assume it is a dCVnet style object
  #   such as dCVnet$final
  if (!"performance" %in% names(x)) {
    stop("not a suitable list for classperformance")
  }
  expected <- c("label", "reference", "probability", "classification")
  if ( any(!expected %in% names(x$performance)) ) {
    stop("required columns for classperformance missing.")
  }
  classperformance.dCVnet(x, ...)
}

#' @export
classperformance.dCVnet <- function(x, as.data.frame = T, ...) {
  if ( as.data.frame ) return(x$performance)
  labs <- as.character(unique(x$performance$label))
  names(labs) <- labs
  R <- lapply(labs,
              function(x) {
                x$performance[x$performance$label == x, ]
              })
  return(structure(R, class = c("classperformance", "list")))
}

#' @export
classperformance.glm <- function(x,
                                 as.data.frame = T,
                                 label = deparse(substitute(x)),
                                 threshold = 0.5, ...) {
  # Return (labelled) prediction dataframe from a glm
  #     given a threshold (default = 0.5):
  lvl <- levels(x$data[, 1])
  classification <- as.numeric(stats::fitted(x) > threshold) + 1
  classification <- factor(lvl[classification], levels = lvl)

  reference <- x$data[, 1]
  probability <- fitted(x)

  R <- data.frame(reference = reference,
                  probability = probability,
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

#' @export
classperformance.glmlist <- function(x, as.data.frame = T, ...) {
  # applies pobj.glm to a list of glms.

  class_list <- c("classperformance", "list")
  class_df <- c("classperformance", "data.frame")

  R <- lapply(1:length(x), function(i) {
    # for a list we force return of a dataframe as we wrap in a list anyway.
    classperformance.glm(x[[i]],
                         as.data.frame = T,
                         label = names(x)[i], ...)
  })
  names(R) <- names(x)
  if ( !as.data.frame ) return(structure(R, class = class_list))
  R <- do.call(rbind, R)
  rownames(R) <- NULL
  return(structure(R, class = class_df))
}


#  ~ using class performance ----------------------------------------------


#' @export
summary.classperformance <- function(object, label = NA, ...) {
  # Function assigns a label if asked to (for multi performance mode.)
  #   In multiperformance mode it uses label to produce columns of data.
  # If label is 'None' then this column is removed.

  if ( "list" %in% class(object) ) {
    # we merge any lists.
    object <- structure(data.frame(do.call(rbind, object)),
                             class = c("classperformance", "data.frame"))
  }

  # Check structure:
  test <- names(object)
  test_cols <- c("reference", "probability", "classification", "label")
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
                           predicted = performance$probability)
    B <- pmax(B, 1 - B)
    B <- data.frame(Measure = "AUROC", Value = B)

    B$label <- A$label <- unique(performance$label)
    return(rbind(A, B))
  }

  .multi_cpsummary <- function(performance) {
    R <- lapply(1:length(unique(performance$label)),
                function(i){
                  rr <- as.character(unique(performance$label)[i])
                  dd <- performance[performance$label == rr, ]
                  R <- summary.classperformance(dd, rr)
                  # Parse back to long.
                  R$label <- NULL # rr
                  names(R)[2] <- rr #"Value"
                  return(R)
                })
    R <- Reduce(function(x, y) merge(x, y, by = "Measure", sort = F), R)
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
  ols <- summary(classperformance(dCVnet_object))

  summary_measures <- c("mean", "sd", "min", "max")
  names(summary_measures) <- summary_measures

  S <- lapply(summary_measures, function(M) {
    apply(ols[, -1], 1, function(i) get(M)(i, na.rm = T))
  } )

  S <- do.call(data.frame, list(Measure = ols[, 1], S))

  S <- data.frame(S, "..." = " - ", ols[, -1], stringsAsFactors = F)

  return(S)
}
