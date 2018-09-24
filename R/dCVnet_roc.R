

# S3 class: rocdata --------
#   acts on classperformance to produce data for a ROC.
#   which is a table of:
#     Sens : Sensitivity (tpr),
#     InvSpec : 1 - Specificity (fpr)
#   for assorted thresholds of the model class thresholds.
# Uses ROCR.


#' extract_rocdata
#'
#' Function converts a classperformance data.frame to the values required for
#'     calculating a ROC curve. It does this by calculating sensitivity and
#'     specificity at incrementing class probability thresholds.
#'
#' Uses \pkg{ROCR} functions: \code{\link[ROCR]{prediction}}, and
#'     \code{\link[ROCR]{performance}}
#'
#' @name extract_rocdata
#' @param classperformance a \code{\link{classperformance}} object.
#'
#' @return a data.frame object of class "rocdata" which can be plotted. Contents:
#'     \itemize{
#'     \item{\code{Sens} : Sensitivity}
#'     \item{\code{InvSpec} : 1 - Specificity}
#'     \item{\code{alpha} : threshold}
#'     \item{\code{run} : label from \code{\link{classperformance}}}
#'     }
#'
#' @seealso \code{\link{plot.rocdata}}
#'
#' @export
extract_rocdata <- function(classperformance) {
  # First ensure it is in list format, not dataframe:
  if ( "data.frame" %in% class(classperformance) ) {
    lvls <- as.character(unique(classperformance$label))
    classperformance <- lapply(lvls,
                          function(lab) {
                            classperformance[classperformance$label == lab, ]
                          })
    names(classperformance) <- lvls
    classperformance <- structure(classperformance, class = c("classperformance", "list"))
  }

  # Utility subfunction convert a ROC performance object to a dataframe.
  .performance_to_data_frame <- function(perf, names) {
    ns <- sapply(perf@y.values, length)
    runs <- rep(names, ns)
    # convert to dataframe for ggplot:
    outer.df <- data.frame(Sens = do.call(c, perf@y.values),
                           InvSpec = do.call(c, perf@x.values),
                           alpha = do.call(c, perf@alpha.values),
                           run = runs)
    return(outer.df)
  }

  outer.pred <- ROCR::prediction(
    predictions = lapply(classperformance, "[[", "probability"),
    labels = lapply(classperformance, "[[", "reference"))
  outer.perf <- ROCR::performance(outer.pred, "tpr", "fpr")

  R <- .performance_to_data_frame(outer.perf, names(classperformance))
  return(structure(R, class = c("rocdata", "data.frame")))
}
