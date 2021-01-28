
# calibration S3 generic ---------------------------------------------

# implements 'weak' (i.e. linear) calibration.


#' calibration
#'
#' calculates 'weak' (i.e. intercept + slope) calibration equivalent to the
#'     values returned by the val.prob function in the rms package
#'     (which accompanies Frank Harrell's Regression Modelling Strategies book).
#'
#' Calibration is not returned via \code{\link{performance}} because of the
#'     computational overhead of model fitting.
#'
#' @name calibration
#'
#' @param x an object containing predicted probabilities and observed outcomes,
#'     from which calibration can be extracted.
#' @param ... arguments to pass on
#'
#' @return calibration intercept and calibration slope
#'
#' @importFrom stats qlogis
#'
#' @export
calibration <- function(x, ...) {
  UseMethod("calibration")
}

# Internal function to run the actual calculation:
calibration_ <- function(predictedprobability,
                         observedoutcome) {
  logit <- qlogis(predictedprobability)
  finite <- is.finite(logit)
  if ( any(!finite) ) {
    # remove from fit:
    logit <- logit[finite]
    observedoutcome <- observedoutcome[finite]
    # warn:
    n <- length(finite)
    i <- sum(!finite)
    warning(paste0(i, "/", n,
                   "values removed due to predictions equal to 0 / 1"))
  }
  # xx is data.frame with reference & prediction (no group)
  m <- glm(formula = observedoutcome ~ logit,
           family = "binomial")
  cc <- coef(m)
  return(c(Intercept = cc[[1]], Slope = cc[[2]]))
}


#' calibration.default
#'
#' @rdname calibration
#' @description Default function behaviour assumes input is a
#'     list/data.frame with required vectors as elements.
#' @param ppid indicator for predicted probability element in x
#'     (column "name" or index)
#' @param ooid indicator for observed outcome element in x
#'     (column "name" or index)
#' @param gid indicator for grouping variable in x
#'     (column "name" or index). Set to NA to force no grouping.
#' @export
calibration.default <- function(x,
                                ppid = "prediction",
                                ooid = "reference",
                                gid = "label",
                                ...) {
  if ( length(unique(c(ppid, ooid, gid))) != 3 ) {
    stop("Predictions and observations must be distinct elements in x")
  }
  pp <- x[[ppid]]
  oo <- x[[ooid]]

  # No implicit recyling
  if ( length(pp) != length(oo) ) {
    stop("numbers of predicted probabilities and
         observed outcomes must match")
  }

  # grouping variable?
  if ( is.na(gid) ) {
    g <- rep(deparse(substitute(x)), times = length(pp))
  } else {
    g <- x[[gid]]
  }

  # iterable for groups:
  glist <- setNames(unique(g), unique(g))

  # calibration for each group:
  result <- vapply(
    glist,
    FUN = function(i) {
      calibration_(predictedprobability = pp[g %in% i],
                   observedoutcome = oo[g %in% i])
    },
    FUN.VALUE = c(Intercept = 0.0,
                  Slope = 1.0),
    USE.NAMES = TRUE)

  class(result) <- c("calcoefs")
  result
}

# Simple print function for calibration result
#' @export
print.calcoefs <- function(x, ...) {
  cat("Calibration coefficents\n")
  print(unclass(x))
  invisible(x)
}

# Simple summary function for calibration result
#' @export
summary.calcoefs <- function(object,
                             FUNS = c(mean = mean,
                                      sd = sd,
                                      min = min,
                                      max = max),
                             ...) {
  if ( NCOL(object) < 2 ) return(object)
  return(vapply(FUNS, function(f) {
    apply(object, 1, f)
  }, FUN.VALUE = 0.0))
}
