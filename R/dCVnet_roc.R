

# S3 class: rocdata --------
#   data for a ROC analysis, acts on performance to produce a table of:
#     Sens : Sensitivity (tpr),
#     InvSpec : one-minus Specificity (fpr)
#   for assorted thresholds of the model class thresholds.
# Uses ROCR.
#


#' extract_rocdata
#'
#' Function reads ROC data from a \code{\link{performance}} object.
#'     Sensitivity and specificity are calculated
#'     at incrementing class probability thresholds.
#'
#' Uses \pkg{ROCR} functions: \code{\link[ROCR]{prediction}}, and
#'     \code{\link[ROCR]{performance}}
#'
#' @name extract_rocdata
#' @param performance a \code{\link{performance}} object.
#' @param invertprob boolean. Should class probabilities be inverted?
#'     If the ROC curve appears under the diagonal then toggle this option.
#' @return a data.frame object of class "rocdata" which can be plotted.
#'     Contents:
#'     \itemize{
#'     \item{\code{Sens} : Sensitivity}
#'     \item{\code{InvSpec} : 1 - Specificity}
#'     \item{\code{alpha} : threshold}
#'     \item{\code{run} : label from \code{\link{performance}}}
#'     }
#'
#' @seealso \code{\link{plot.rocdata}}
#'
#' @export
extract_rocdata <- function(performance,
                            invertprob = FALSE) {
  # First ensure it is in list format, not dataframe:
  if ( "data.frame" %in% class(performance) ) {
    lvls <- as.character(unique(performance$label))

    performance <- lapply(
      lvls,
      function(lab) {
        performance[performance$label == lab, ]
      })
    names(performance) <- lvls

    performance <- structure(performance,
                                  class = c("performance", "list"))
  }

  # Utility subfunction convert a ROC performance object to a dataframe.
  .performance_to_data_frame <- function(perf, names) {
    ns <- vapply(perf@y.values, length, c(1.0))
    runs <- rep(names, ns)
    # convert to dataframe for ggplot:
    outer.df <- data.frame(Sens = do.call(c, perf@y.values),
                           InvSpec = do.call(c, perf@x.values),
                           alpha = do.call(c, perf@alpha.values),
                           run = runs,
                           stringsAsFactors = FALSE)
    return(outer.df)
  }

  .extract_pred <- function(x, invert) {
    if ( invert ) {
      return(1 - x$prediction)
    } else {
      return(x$prediction)
    }
  }

  outer.pred <- ROCR::prediction(
    # Invert class probabilities if needed:
    predictions = lapply(performance,
                         .extract_pred,
                         invert = invertprob),
    labels = lapply(performance, "[[", "reference"),
    label.ordering = levels(performance[[1]]$reference))
  outer.perf <- ROCR::performance(outer.pred, "tpr", "fpr")

  R <- .performance_to_data_frame(outer.perf, names(performance))
  return(structure(R, class = c("rocdata", "data.frame")))
}


#' average_rocdata
#'
#' Function calculates the average of a set of logistic regression ROC curves.
#'     For a range of thresholds between 0 and 1 the sensitivity and
#'     specificity are extracted and mean averaged over the set of curves.
#'     This is intended to be used with k-fold cross-validation data.
#'
#'     Warning: The averaging of roc-curves is difficult
#'         (e.g. \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4112395/}{
#'         Chen & Samuelson Br J. Radiol Aug 2014})
#'         , particularly while preserving the AUROC such that the AUROC
#'         for the average ROC-curve is equal to the average AUROC of the
#'         component ROC-curves.
#'
#'         This function is for a rough display of average performance
#'         over bootstraps or repeats of k-fold cross-validation.
#'
#'         Minimal testing suggests that there should be agreement in the
#'         AUROCs to the third decimal place.
#'
#' @name average_rocdata
#' @param rocdata a \code{\link{extract_rocdata}} object.
#' @param n the number of thresholds in \\[0,1\\] to evaluate.
#'
#' @return a data.frame object of class "rocdata" which can be plotted.
#'     Contents:
#'     \itemize{
#'     \item{\code{Sens} : Sensitivity}
#'     \item{\code{InvSpec} : 1 - Specificity}
#'     \item{\code{alpha} : threshold}
#'     \item{\code{run} : a label indicating this is averaged}
#'     }
#'
#' @seealso \code{\link{plot.rocdata}}, \code{\link{extract_rocdata}}
#'
#' @export
average_rocdata <- function(rocdata,
                            n = 1000) {
  # rocdata is output of extract_rocdata.

  # 1 roc curve per level of rocdata$run:
  runs <- as.character(unique(rocdata$run))
  # Split data to a list of rocdatas:
  ds <- lapply(runs, function(x) rocdata[rocdata$run == x, ] )

  # Calculate sens and spec values
  #     at a set of standard thresholds [0,1] for each curve.

  # set the fixed alphas (Inf needs special treatment because of findInterval)
  alphas <- c(seq(0, 1, length.out = n - 1), Inf)
  alphalabs <- seq_along(alphas)

  res <- lapply(ds, function(d) {
    d <- d[order(d$alpha), ]
    dinf <- d[is.infinite(d$alpha), ] # set aside Inf

    d <- d[findInterval(alphas, vec = d$alpha,
                        all.inside = TRUE, left.open = TRUE), ]

    d[is.infinite(alphas), ] <- dinf  # add inf back in.

    d$alab <- alphalabs # this label is used as the grouping variable
    return(d)
  } )
  # merge results:
  res <- as.data.frame(data.table::rbindlist(res), stringsAsFactors = FALSE)

  # average results for threshold values:
  av <- aggregate(res[, !names(res) %in% c("run", "alab")],
                  by = list(alab = res$alab),
                  FUN = "mean")
  # add a run label
  av$run <- "Average"

  av <- av[order(av$alpha, decreasing = TRUE), -1]
  av <- unique(av)
  rownames(av) <- NULL
  # Use this to return a combined dataset: rbind(rocdata, av)
  return(structure(av, class = c("rocdata", "data.frame")))
}
