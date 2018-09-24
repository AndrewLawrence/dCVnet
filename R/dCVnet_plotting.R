
# Plotting Functions ------------------------------------------------------


#  ~ for S3 classes -------------------------------------------------------

#' @export
plot.repeated.cv.glmnet <- function(x, ...) {
  x$lambda <- log10(x$lambda)
  p <- ggplot2::ggplot(x,
                       ggplot2::aes_string(y = "cvm", x = "lambda")) +
    ggplot2::geom_line(data = x, colour = "grey", size = 1.2) +
    ggplot2::geom_point(data = x[x$lambda.min, ],
                        ggplot2::aes_string(y = "cvm", x = "lambda"),
                        shape = 25, size = 2, colour = "red", fill = "black",
                        inherit.aes = F) +
    ggplot2::ylab(paste("Metric:", attr(x, "type.measure"))) +
    ggplot2::xlab("lambda (log10)") +
    ggplot2::theme_light()
  print(p)

  invisible(x)
}

#' @export
plot.multialpha.repeated.cv.glmnet <- function(x, ...) {
  x <- x$inner_results # only need the inner_results dataframe.

  p <- ggplot2::ggplot(x,
                       ggplot2::aes_string(y = "cvm", x = "log10(lambda)",
                                    colour = "alpha", fill = "alpha")) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = x[x$lambda.min, ],
                        colour = "black", shape = 25) +
    ggplot2::geom_point(data = x[x$best, ],
                        colour = "black", shape = 1, size = 5,
                        show.legend = F) +
    ggplot2::ylab(paste("Metric:", attr(x, "type.measure"))) +
    ggplot2::xlab("Lambda path (log10)") +
    ggplot2::theme_light()
  print(p)

  invisible(x)
}

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


tuning_plot_dCVnet <- function(object) {
  # Plotting function to show outer fold variability in the tuning curves
  #   at different alphas.

  # Merge datasets adding foldids:
  df <- lapply(1:length(object$tuning),
               function(x) {
                 R <- object$tuning[[x]]$tuning$inner_results
                 R$outfold <- names(object$tuning)[x]
                 return(R)
               })

  df <- do.call(rbind, df)
  rownames(df) <- NULL

  df$Fold <- sapply(strsplit(df$outfold, split = "\\."), "[", 1)
  df$Rep <- sapply(strsplit(df$outfold, split = "\\."), "[", 2)

  p <- ggplot2::ggplot(df,
                       ggplot2::aes_string(y = "cvm", x = "log10(lambda)",
                                    colour = "alpha", fill = "alpha",
                                    label = "gsub(\"OutFold\", \"\", Fold)",
                                    group = "paste(alpha, outfold)")) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = df[df$lambda.min, ],
                        colour = "black", size = 3, shape = 24) +
    ggplot2::geom_text(data = df[df$lambda.min, ],
                       colour = "black", size = 3) +
    ggplot2::ylab(paste("Metric:", attr(df, "type.measure"))) +
    ggplot2::xlab("lambda path (log10)") +
    ggplot2::facet_wrap(~Rep) +
    ggplot2::theme_light()
  print(p)
  invisible(df)
}


# Visualise:
#' plot_outerloop_coefs
#'
#' Plot showing standardised betas for model coefficients.
#'     coefficients are first mean-averaged over the k-folds,
#'     then displayed per-repetition.
#'
#' @param object a dCVnet object
#' @importFrom scales muted
#' @export
plot_outerloop_coefs <- function(object) {
  df <- coef.dCVnet(object, type = "rep")
  df$StdBeta <-  df$Coef

  p <-  ggplot2::ggplot(df,
                        ggplot2::aes_string(y = "StdBeta",
                                            x = "Predictor",
                                            colour = "StdBeta")) +
    ggplot2::geom_hline(yintercept = 0.0) +
    ggplot2::geom_point(alpha = 1) +
    ggplot2::scale_colour_gradient2(low = muted("blue"),
                                    mid = "grey",
                                    high = muted("red"),
                                    midpoint = 0) +
    ggplot2::stat_summary(fun.data = "mean_se") +
    ggplot2::ylab("Standardised Beta") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1.0,
                                                       vjust = 0.5)) +
    ggplot2::theme_light()
  print(p)

  invisible(df)
}
