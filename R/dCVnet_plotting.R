
# Plotting Functions ------------------------------------------------------


#  ~ for S3 classes -------------------------------------------------------

#' plot.multialpha.repeated.cv.glmnet
#'
#' Plot averaged and best innerloop model performance highlighting
#'     the final selected alpha for a
#'     \code{\link{multialpha.repeated.cv.glmnet}}
#'
#' @param x a \code{\link{multialpha.repeated.cv.glmnet}} object.
#' @param xvar what to plot on x-axis. Either the log10 lambda value, or the
#'     proportion of the lambda sequence (for that alpha).
#' @param errorbars boolean. Add errorbars to the plot?
#' @param ... NULL
#'
#' @inheritParams tuning_plot_dCVnet
#'
#' @export
plot.multialpha.repeated.cv.glmnet <- function(x,
                                               xvar = c("lambda", "s"),
                                               errorbars = FALSE,
                                               plot = TRUE,
                                               ...) {
  xvar <- match.arg(xvar)
  # store attributes:
  type.lambda <- attr(x, "type.lambda")
  type.measure <- attr(x, "type.measure")

  x <- x$results
  # proc data:
  x$alpha <- as.factor(x$alpha)
  if ( identical(xvar, "s") ) { # nolint
    x$s <- as.numeric(gsub("s", "", x$s))

    smax <- aggregate(x$s, by = list(alpha = x$alpha), FUN = max)

    x$s <- x$s / smax$x[match(x$alpha, smax$alpha)]
  }

  p <- ggplot2::ggplot(x,
                       ggplot2::aes_string(y = "cvm",
                                           ymin = "cvlo",
                                           ymax = "cvup",
                                           x = xvar,
                                           colour = "alpha",
                                           fill = "alpha")) +
    ggplot2::geom_line()

  # error-bars?
  if ( errorbars ) {
    p <- p + ggplot2::geom_linerange()
  }

  p <- p +
    ggplot2::geom_point(shape = 6) +
    ggplot2::geom_point(data = x[x[[type.lambda]], ],
                        colour = "black", shape = 25) +
    ggplot2::geom_point(data = x[x$best, ],
                        colour = "black", shape = 1, size = 5,
                        show.legend = FALSE) +
    ggplot2::ylab(paste("Metric:", type.measure)) +
    ggplot2::theme_light()

  # x-axis options:
  if ( identical(xvar, "lambda") ) {
    p <- p +
      ggplot2::scale_x_log10() +
      ggplot2::xlab(paste0("Lambda (log10)\nSelection: ", type.lambda))
  } else {
    p <- p +
      ggplot2::xlab(paste0("Lambda path fraction\nSelection: ", type.lambda))
  }
  if ( plot ) print(p)

  invisible(list(plot = p,
                 data = x))
}

#' plot.dCVnet
#'
#' Produces either a plot of the tuning parameters and inner loop performance
#'     per fold of the outer loop, or the final roc plots of the outer loop.
#'     \itemize{
#'     \item{\code{tuning} - selected tuning parameters over repeated folds
#'         of the outerloop}
#'     \item{\code{roc} - Sensitivity vs. (1 - Specificity) plot.}
#'     }
#'
#' @name plot.dCVnet
#' @param x a \code{\link{dCVnet}} object
#' @param type one of
#'     \itemize{
#'     \item{\code{"tuning"} - hyperparameter tuning plots
#'          (\code{\link{tuning_plot_dCVnet}})}
#'     \item{\code{"roc"} - a roc plot (\code{\link{plot.rocdata}})}
#'     }
#' @param ... additional arguments passed to plot functions above.
#' @export
plot.dCVnet <- function(x, type = "tuning", ...) {
  # options:
  #   "tuning"
  #   "ROC"
  type <- match.arg(type, choices = c("tuning", "roc"), several.ok = FALSE)

  f <- family(x)
  if ( type == "roc" && !(f %in% c("binomial"))) {
    stop("roc plots supported for binomial only")
  }

  switch(type,
         tuning = return(tuning_plot_dCVnet(x, ...)),
         roc = return(plot(extract_rocdata(performance(x)), ...))
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
#' @param legend logical. Display a legend?
#' @param alphalabel should certain alpha values (probability thresholds)
#'     on the curve be highlighted with symbols indicating threshold?
#' @param guide_labels a named list of 3 labels used in the legend:
#'     \itemize{
#'     \item{group = label for the Grouping factor}
#'     \item{threshold = label for the Threshold factor}
#'     \item{refline = label for the Reference line}
#'     }
#' @param ... additional arguments (unused)
#' @return a ROC plot, as above.
#'
#' @export
plot.rocdata <- function(x,
                         legend = TRUE,
                         alphalabel = c(0.25, 0.5, 0.75),
                         guide_labels = c(group = "Model",
                                          threshold = expression(P[Threshold]),
                                          refline = "Chance\nPerformance"),
                         ...) {
  guide_labels_expect <- c("group", "threshold", "refline")
  if ( ! all( guide_labels_expect %in% names(guide_labels)) ) {
    warning(paste("guide_labels must contain named elements:",
                  "\t\t'group', 'threshold' and 'refline'",
                  " Using defaults instead", sep = "\n"))
    guide_labels <- c(group = "Model",
                      threshold = expression(P[Threshold]),
                      refline = "Chance\nPerformance")
  }
  # were thresholds supplied?
  hasalphas <- !any(is.na(alphalabel)) && !identical(FALSE, alphalabel)

  .closest <- function(vals, x) {
    locs <- vapply(vals, function(v) which.min(abs(x - v)), FUN.VALUE = 1L)
    r <- rep(NA_character_, NROW(x))
    r[locs] <- names(locs)
    return(r)
  }

  if ( hasalphas ) {
    alphalabel <- as.numeric(alphalabel)
    stopifnot(min(alphalabel) >= 0.0,
              max(alphalabel) <= 1.0)

    if ( is.null(names(alphalabel)) ) {
      alphalabel <- setNames(alphalabel, prettyNum(alphalabel))
    }
    # initialise in x:
    x$PThreshold <- NA_character_
    # get closest for each group (called run):
    # note: lappy used for side effect.
    lapply(
      setNames(unique(x$run),
               unique(x$run)),
      function(g) {
        x$PThreshold[x$run == g] <<- .closest(alphalabel, x$alpha[x$run == g])
      })
    x$PThreshold <- factor(x$PThreshold, levels = names(alphalabel))
  }

  refline <- data.frame(Sens = 0,
                        InvSpec = 0,
                        SensEnd = 1,
                        InvSpecEnd = 1,
                        lty = guide_labels$refline,
                        stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(x,
                       ggplot2::aes_string(y = "Sens",
                                           x = "InvSpec")) +
    ggplot2::geom_segment(
      ggplot2::aes_(
        x = ~InvSpec,
        y = ~Sens,
        xend = ~InvSpecEnd,
        yend = ~SensEnd,
        lty = ~lty),
      show.legend = legend,
      data = refline,
      inherit.aes = FALSE) +
    ggplot2::geom_line(ggplot2::aes_string(colour = "run"),
                       show.legend = legend) +
    ggplot2::xlab("False Positive Rate") +
    ggplot2::ylab("True Positive Rate") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1), minor_breaks = NULL) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), minor_breaks = NULL) +
    ggplot2::scale_size_manual(values = c(0.7)) +
    ggplot2::coord_equal() +
    ggplot2::theme_light()

  if ( legend ) {
    p <- p + ggplot2::guides(
      colour = ggplot2::guide_legend(title = guide_labels$group, order = 1),
      lty = ggplot2::guide_legend(title = NULL,
                                  order = 3,
                                  override.aes = list(
                                    linetype = 1,
                                    shape = 32,
                                    alpha = 1)))
  }

  if ( hasalphas ) {
    p <- p +
      ggplot2::geom_point(data = x[!is.na(x$PThreshold), ],
                          mapping = ggplot2::aes_string(shape = "PThreshold",
                                                        colour = "run"),
                          show.legend = legend)

    if ( legend ) {
      p <- p + ggplot2::guides(
        shape = ggplot2::guide_legend(title = guide_labels$threshold,
                                      order = 2,
                                      override.aes = list(linetype = NA))
      )
    }
  }
  print(p)
  return(invisible(list(plot = p,
                        data = x)))
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
#' @param plot (bool) produce the plot? or just return the plot and data
#' @return a list containing the plot and a data.frame used to plot
#'     (full data.frame is returned, i.e. ignores n.random)
#'
#' @export
tuning_plot_dCVnet <- function(object, n.random = 0, plot = TRUE) {
  # Plotting function to show outer fold variability in the tuning curves
  #   at different alphas.

  # Merge datasets adding foldids:
  df <- lapply(seq_along(object$tuning),
               function(x) {
                 R <- object$tuning[[x]]$tuning$results
                 R$outfold <- names(object$tuning)[x]
                 return(R)
               })

  df <- as.data.frame(data.table::rbindlist(df), stringsAsFactors = FALSE)
  rownames(df) <- NULL

  df$Fold <- vapply(X = strsplit(df$outfold, split = "\\."),
                    FUN = "[",
                    FUN.VALUE = c(""),
                    1)
  df$Rep <- vapply(X = strsplit(df$outfold, split = "\\."),
                   FUN = "[",
                   FUN.VALUE = c(""),
                   2)

  plotReps <- unique(df$Rep)
  if ( n.random > 0 ) plotReps <- sample(plotReps, size = n.random)
  pdf <- df[df$Rep %in% plotReps, ]
  # should we treat alpha as continuous or discrete?
  if ( length(unique(pdf$alpha)) < 7 ) pdf$alpha <- as.factor(pdf$alpha)

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

  if ( plot ) print(p)
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
#' @param ordered sort predictors by size?
#' @param abs plot absolute values?
#' @param intercept include the value of the intercept coefficient in the plot?
#' @param final add the final model coefficients as an overlay?
#'     (note: this is not in the return data, but can be accessed with
#'      coef(object$final$model))
#' @param final_col colour for final model coefficients
#' @param final_shape shape for final model coefficients
#' @importFrom scales muted
#' @export
plot_outerloop_coefs <- function(object,
                                 type = "rep",
                                 ordered = FALSE,
                                 abs = FALSE,
                                 intercept = FALSE,
                                 final = TRUE,
                                 final_col = "red",
                                 final_shape = 24) {
  df <- coef.dCVnet(object, type = type)

  if ( ! intercept ) {
    df <- df[df$Predictor != "(Intercept)", ]
  }

  if ( abs ) {
    df$stdbeta <- abs(df$Coef)
    ylabel <- "|Standardised Beta|"
  } else {
    df$stdbeta <- df$Coef
    ylabel <- "Standardised Beta"
  }

  if ( ordered ) {
    df$Predictor <- forcats::fct_reorder(df$Predictor, df$stdbeta)
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

  if ( final ) {
    fcoef <- coef(object$final$model)
    if ( abs ) {
      fcoef <- abs(fcoef)
    }
    fcoef.df <- data.frame(Predictor = factor(rownames(fcoef),
                                              levels = levels(df$Predictor)),
                           stdbeta = as.vector(fcoef),
                           fold = "Final",
                           rep = "Final",
                           stringsAsFactors = FALSE)
    if ( ! intercept ) {
      fcoef.df <- fcoef.df[fcoef.df$Predictor != "(Intercept)", ]
    }

    p <- p + ggplot2::geom_point(data = fcoef.df,
                                 colour = final_col,
                                 shape = final_shape)
  }

  print(p)

  invisible(list(plot = p,
                 data = df))
}
