

# S3 class: rocdata --------
#   acts on classperformance to produce data for a ROC.
#   which is a table of:
#     Sensitivity (tpr),
#     1-Specificity (fpr)
#   for assorted thresholds of the model class thresholds.
# Uses ROCR.


rocdata <- function(x, ...) {
  UseMethod("rocdata", x)
}

rocdata.classperformance <- function(performance) {
  # First ensure it is in list format, not dataframe:
  if ( "data.frame" %in% class(performance) ) {
    lvls <- as.character(unique(performance$label))
    performance <- lapply(lvls,
                          function(lab) {
                            performance[performance$label == lab,]
                          })
    names(performance) <- lvls
    performance <- structure(performance, class = c("classperformance", "list"))
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
    predictions = lapply(performance, "[[", "probability"),
    labels = lapply(performance, "[[", "reference"))
  outer.perf <- ROCR::performance(outer.pred, "tpr", "fpr")

  R <- .performance_to_data_frame(outer.perf, names(performance))
  return(structure(R, class = c("rocdata", "data.frame")))
}



