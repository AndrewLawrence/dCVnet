

# classperformance S3 generic ---------------------------------------------


# a dCVnet S3 generic: classperformance
#   returns a merged dataframe or a list of data frames
#   where the dataframe contains raw performance information.


#  ~ making classperformance ----------------------------------------------


classperformance <- function(x, as.data.frame = T, ...) {
  UseMethod("classperformance", x)
}


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


classperformance.dCVnet <- function(object, as.data.frame = T) {
  if ( as.data.frame ) { return(object$performance) }
  labs <- as.character(unique(object$performance$label))
  names(labs) <- labs
  R <- lapply(labs,
              function(x) {
                object$performance[object$performance$label == x,]
              })
  return(structure(R, class = c("classperformance", "list")))
}


classperformance.glm <- function(glm,
                                 as.data.frame = T,
                                 label = deparse(substitute(glm)),
                                 threshold = 0.5) {
  # Return (labelled) prediction dataframe from a glm
  #     given a threshold (default = 0.5):
  lvl <- levels(glm$data[,1])
  classification <- as.numeric(fitted(glm) > threshold) + 1
  classification <- factor(lvl[classification], levels = lvl)

  reference <- glm$data[,1]
  probability <- fitted(glm)

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


classperformance.glmlist <- function(glmlist, as.data.frame = T, ...) {
  # applies pobj.glm to a list of glms.
  R <- lapply(1:length(glmlist), function(i) {
    # for a list we force return of a dataframe as we wrap in a list anyway.
    classperformance.glm(glmlist[[i]],
                         as.data.frame = T,
                         label = names(glmlist)[i], ...)
  })
  names(R) <- names(glmlist)
  if ( !as.data.frame ) { return(structure(R, class = c("classperformance", "list"))) }
  R <- do.call(rbind, R)
  rownames(R) <- NULL
  return(structure(R, class = c("classperformance", "data.frame")))
}


#  ~ using class performance ----------------------------------------------



summary.classperformance <- function(performance, label = NA) {
  # Function assigns a label if asked to (for multi performance mode.)
  #   In multiperformance mode it uses label to produce columns of data.
  # If label is 'None' then this column is removed.

  if ( "list" %in% class(performance) ) {
    # we merge any lists.
    performance <- structure(data.frame(do.call(rbind, performance)),
                             class = c("classperformance", "data.frame"))
  }

  # Check structure:
  test <- names(performance)
  test <- any(!(c("reference", "probability", "classification", "label") %in% test))
  if ( test ) { cat(names(performance)); stop("Check input") }

  if ( !is.na(label) ) { performance$label <- label }

  # Two methods:
  .singleclassperformance_summary <- function(performance) {
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

  .multiclassperformance_summary <- function(performance) {
    R <- lapply(1:length(unique(performance$label)),
                function(i){
                  rr <- as.character(unique(performance$label)[i])
                  dd <- performance[performance$label == rr,]
                  R <- summary.classperformance(dd, rr)
                  # Parse back to long.
                  R$label <- NULL # rr
                  names(R)[2] <- rr #"Value"
                  return(R)
                })
    R <- Reduce(function(x, y) { merge(x, y, by = "Measure", sort = F) }, R)
    return(R)
  }
  # If we lack a label col assign it and return
  #   (for multiclassperf functionality)
  # If there is a single label return single, otherwise return mulit.
  if ( is.null(performance$label) ) {
    R <- .singleclassperformance_summary(performance)
  } else {
    if ( length(unique(performance$label)) == 1 ) {
      R <- .singleclassperformance_summary(performance)
    } else {
      R <- .multiclassperformance_summary(performance)
    }
  }
  # Option : remove label if it is there.
  if ( label %in% "None" ) { R$label <- NULL }
  return(R)
}


report_classperformance_summary <- function(object) {
  ols <- summary(classperformance(object))

  summary_measures <- c("mean", "sd", "min", "max")
  names(summary_measures) <- summary_measures

  S <- lapply(summary_measures, function(M) {
    apply(ols[,-1], 1, function(i) get(M)(i, na.rm = T))
  } )

  S <- do.call(data.frame, list(Measure = ols[,1], S))

  S <- data.frame(S, "..." = " - ", ols[,-1], stringsAsFactors = F)

  return(S)
}

