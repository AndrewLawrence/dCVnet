
# Plotting Functions ------------------------------------------------------


#  ~ for S3 classes -------------------------------------------------------

#' plot.repeated.cv.glmnet
#'
#' Plot averaged and best innerloop model performance for a (single alpha)
#'     \code{\link{repeated.cv.glmnet}}
#'
#' @param x a \code{\link{repeated.cv.glmnet}} object.
#' @param ... NULL
#' @export
plot.repeated.cv.glmnet <- function(x, ...) {
  p <- ggplot2::ggplot(x,
                       ggplot2::aes_string(y = "cvm", x = "lambda")) +
    ggplot2::geom_line(data = x, colour = "grey", size = 1.2) +
    ggplot2::geom_point(data = x[x$lambda.min, ],
                        ggplot2::aes_string(y = "cvm", x = "lambda"),
                        shape = 25, size = 2, colour = "red", fill = "black",
                        inherit.aes = F) +
    ggplot2::scale_x_log10() +
    ggplot2::ylab(paste("Metric:", attr(x, "type.measure"))) +
    ggplot2::xlab("lambda (log10)") +
    ggplot2::theme_light()
  print(p)

  invisible(p)
}

#' plot.repeated.cv.glmnet
#'
#' Plot averaged and best innerloop model performance highlighting
#'     the final selected alpha for a
#'     \code{\link{multialpha.repeated.cv.glmnet}}
#'
#' @param x a \code{\link{multialpha.repeated.cv.glmnet}} object.
#' @param ... NULL
#' @export
plot.multialpha.repeated.cv.glmnet <- function(x, ...) {
  x <- x$inner_results # only need the inner_results dataframe.

  p <- ggplot2::ggplot(x,
                       ggplot2::aes_string(y = "cvm", x = "lambda",
                                           colour = "alpha", fill = "alpha")) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = x[x$lambda.min, ],
                        colour = "black", shape = 25) +
    ggplot2::geom_point(data = x[x$best, ],
                        colour = "black", shape = 1, size = 5,
                        show.legend = F) +
    ggplot2::scale_x_log10() +
    ggplot2::ylab(paste("Metric:", attr(x, "type.measure"))) +
    ggplot2::xlab("Lambda path (log10)") +
    ggplot2::theme_light()
  print(p)

  invisible(list(plot = p,
                 data = x))
}

#' plot.dCVnet
#'
#' Produces either a plot of the tuning parameters and inner loop performance
#'     per fold of the outer loop, or the final roc plots of the outer loop.
#'     \itemize{
#'     \item{\code{tuning} - selected tuning parameters over repeated folds of the
#'           outerloop}
#'     \item{\code{roc} - Sensitivity vs. (1 - Specificity) plot.}
#'     }
#'
#' @name plot.dCVnet
#' @param x a \code{\link{dCVnet}} object
#' @param type one of
#'     \itemize{
#'     \item{\code{"tuning"} - outer tuning coefficients (\code{\link{plot_outerloop_coefs}})}
#'     \item{\code{"roc"} - a roc plot (\code{\link{plot.rocdata}})}
#'     }
#' @param ... additional arguments passed to plot functions above.
#' @export
plot.dCVnet <- function(x, type = "tuning", ...) {
  # options:
  #   "tuning"
  #   "ROC"
  in_type <- type
  type_opts <- c("tuning", "roc")
  type <- pmatch(tolower(type), type_opts)
  if (any(is.na(type))) {
    stop(paste("type: ",
               in_type,
               "must be one of:",
               paste0(type_opts, collapse = ", ")))
  }
  type <- type_opts[type]
  switch(type,
         tuning = return(tuning_plot_dCVnet(x, ...)),
         roc = return(plot(extract_rocdata(classperformance(x)), ...))
  )
}

#' plot.rocdata
#'
#' acts on a \code{\link[=extract_rocdata]{rocdata}} object
#'     to provide a ROC plot with the following:
#'     \itemize{
#'     \item{y-axis : Sensitivity (i.e. True Positive Rate)}
#'     \item{x-axis : 1 - Specificity (i.e. False Positive Rate)}
#'     \item{curve : shows c.d.f. as classification threshold varies}
#'     }
#'
#' requires \pkg{ggplot2}.
#'
#' @name plot.rocdata
#' @param x \code{\link[=extract_rocdata]{rocdata}} object
#' @param legend logical. Display legend?
#' @param ... additional arguments
#' @return a ROC plot, as above.
#'
#' @export
plot.rocdata <- function(x, legend = F, ...) {

  p <- ggplot2::ggplot(x, ggplot2::aes_string(y = "Sens",
                                              x = "InvSpec",
                                              group = "run",
                                              colour = "run")) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = "black") +
    ggplot2::geom_line(show.legend = legend) +
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1), minor_breaks = NULL) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), minor_breaks = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_light()

  return(p)
}


#  ~ non-S3 plotters ------------------------------------------------------


#' tuning_plot_dCVnet
#'
#' Display the tuning curves for the repeated k-folds of the outerloop of a
#' dCVnet.
#'
#' requires \pkg{ggplot2}.
#'
#' @name tuning_plot_dCVnet
#' @param object a \code{\link{dCVnet}} object
#' @param n.random select a random sample of k-fold reps to display.
#'      0 = display all.
#' @return a data.frame containing the full dataset used to plot
#'     (ignores n.random)
#'
#' @export
tuning_plot_dCVnet <- function(object, n.random = 0) {
  # Plotting function to show outer fold variability in the tuning curves
  #   at different alphas.

  # Merge datasets adding foldids:
  df <- lapply(seq_along(object$tuning),
               function(x) {
                 R <- object$tuning[[x]]$tuning$inner_results
                 R$outfold <- names(object$tuning)[x]
                 return(R)
               })

  df <- do.call(rbind, df)
  rownames(df) <- NULL

  df$Fold <- sapply(strsplit(df$outfold, split = "\\."), "[", 1)
  df$Rep <- sapply(strsplit(df$outfold, split = "\\."), "[", 2)

  plotReps <- unique(df$Rep)
  if ( n.random > 0 ) plotReps <- sample(plotReps, size = n.random)
  pdf <- df[df$Rep %in% plotReps, ]

  p <- ggplot2::ggplot(pdf,
                       ggplot2::aes_string(
                         y = "cvm",
                         x = "log10(lambda)",
                         colour = "alpha", fill = "alpha",
                         label = "gsub(\"OutFold\", \"\", Fold)",
                         group = "paste(alpha, outfold)")) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = pdf[pdf$lambda.min, ],
                        colour = "black", size = 3, shape = 24) +
    ggplot2::geom_text(data = pdf[pdf$lambda.min, ],
                       colour = "black", size = 3) +
    ggplot2::ylab(paste("Metric:", attr(pdf, "type.measure"))) +
    ggplot2::xlab("lambda path (log10)") +
    ggplot2::facet_wrap(~Rep) +
    ggplot2::theme_light()
  print(p)
  invisible(list(plot = p, data = df))
}


#' plot_outerloop_coefs
#'
#' Plot showing standardised betas for outer-loop model coefficients.
#'      Because each fold/repetition of the outerloop can have
#'      completely different amounts and types of regularisation
#'      this should be interpreted with caution.
#'
#' It is circular to use these plots to select a subset of variables.
#'      Because these coefficients are based on the complete dataset,
#'      this can produce optimism.
#'
#' @param object a \code{\link{dCVnet}} object
#' @param type How to display coefficients.
#'                 passed to \code{\link{coef.dCVnet}} .
#'     \itemize{
#'     \item{\code{"all"} - boxplot of coefficients for each rep/fold.}
#'     \item{\code{"rep"} - boxplot of mean coefficients for each rep
#'         (mean average over folds).}
#'     \item{\code{"mean"} - dotplot of the mean of \code{"rep"}.}
#'     \item{\code{"median"} - dotplot of the median of \code{"rep"}.}
#'     }
#' @param abs plot absolute values?
#' @importFrom scales muted
#' @export
plot_outerloop_coefs <- function(object, type = "rep", abs = FALSE) {
  df <- coef.dCVnet(object, type = type)

  if ( abs ) {
    df$stdbeta <- abs(df$Coef)
    ylabel <- "|Standardised Beta|"
  } else {
    df$stdbeta <- df$Coef
    ylabel <- "Standardised Beta"
  }

  p <-  ggplot2::ggplot(df,
                        ggplot2::aes_string(y = "stdbeta",
                                            x = "Predictor")) +
    ggplot2::geom_hline(yintercept = 0.0) +
    # ggplot2::scale_fill_gradient2(low = scales::muted("blue"),
    #                                 mid = "grey",
    #                                 high = scales::muted("red"),
    #                                 midpoint = 0) +
    #ggplot2::stat_summary(fun.data = "mean_se") +
    ggplot2::ylab(ylabel) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1.0,
                                                       vjust = 0.5))

  if ( type %in% c("rep", "all") ) {
    p <- p + ggplot2::geom_boxplot()
  } else {
    p <- p + ggplot2::geom_point()
  }

  print(p)

  invisible(list(plot = p,
                 data = df))
}
