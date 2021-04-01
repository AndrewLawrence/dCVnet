
#' performance_with_permuted_labels
#'
#' Reruns a dCVnet model with permuted y-values. This breaks any link between
#'     outcome (y) and features. These results can be used in a permutation test
#'     to verify that dCVnet's double (nested) cross-validation is not
#'     leaking information producing an over-estimate of performance.
#'
#' @section Warning:
#' This code will repeatedly run dCVnet (an already time-consuming program).
#' It should be expected to be *extremely* time-consuming.
#'
#' Further, in most cases evaluating permuted performance is not necessary.
#' Assessing performance with permuted labels is only intended to confirm
#' that the implementation of cross-validation in dCVnet is not leaky.
#'
#' @name performance_with_permuted_labels
#'
#' @param x a dCVnet model.
#' @param n_times number of times to permute labels and obtain
#'     performance measures
#' @return a named list containing the observed (original) performance measures
#'     extracted from x, and a table of equivalent measures obtained under
#'     permutation of the outcome variable.
#' @examples
#' \dontrun{
#'
#' siris <- droplevels(subset(iris, iris$Species != "versicolor"))
#' siris[,1:4] <- scale(siris[,1:4])
#' set.seed(1)
#'
#' model <- dCVnet(y = siris$Species,
#'                 f = ~ Sepal.Length + Sepal.Width +
#'                       Petal.Length + Petal.Width,
#'                 data = siris, nrep_outer = 3, nrep_inner = 3,
#'                 alphalist = c(0.5),
#'                 opt.lambda.type = "1se")
#' nullresult <- performance_with_permuted_labels(model, n_times = 3)
#'
#' # consider the AUROC for the null:
#' range(unlist(nullresult$permuted["AUROC",]))
#' # [1] 0.3627333 0.5235333
#'
#' # vs. the observed value:
#' nullresult$observed["AUROC"]
#' # AUROC
#' #     1
#'
#' }
#'
#' @export
performance_with_permuted_labels <- function(x, n_times = 5) {

  if ( ! x$input$callenv$family %in% c("binomial") ) {
    # currently the function works only with vector y.
    stop("family not supported")
  }

  # function to re-run a dCVnet object with the original data and arguments:
  .rerun_dCVnet <- function(x, shuffley = TRUE) {
    bits <- x$input$callenv
    if ( shuffley ) {
      bits$y <- sample(bits$y)
    }
    do.call("dCVnet", args = bits[-1])
  }

  # extract the parts of performance we are interested in.
  .get_performance <- function(x) {
    p <- performance(x)
    # fix for models estimated in older versions of dCVnet:
    if ( ! "performance" %in% p ) class(p) <- c("performance", "data.frame")
    s <- summary(p)
    s <- setNames(rowMeans(s[, -1]), s[, 1])
    cal <- suppressWarnings(calibration(p))
    cal <- setNames(rowMeans(cal), rownames(cal))
    return(c(s, cal))
  }
  val_labs <- paste0("null", seq_len(n_times))

  res <- parallel::mclapply(setNames(val_labs, val_labs),
                            mc.cores = getOption("mc.cores", 1L),
                            function(i) {
                              .get_performance(.rerun_dCVnet(x,
                                                             shuffley = TRUE))
                            } )
  return(list(observed = .get_performance(x),
              permuted = as.data.frame(res)))
}
