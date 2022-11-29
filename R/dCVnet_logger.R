
# This part of the code is dependent on "openxlsx" package, but this is not a
#   dependency of the package.


#' log_results_to_excel
#'
#' write the results of a \code{\link{dCVnet}} model to a multi-worksheet
#'     xlsx format excel file.
#'
#' @param object a \code{\link{dCVnet}} object.
#' @param file path to a writable excel file (i.e. must not be open/locked).
#' @param referencemodel either set to \code{TRUE} to calculate
#'      reference models, or pass a pre-calculated \code{\link{refunreg}}
#'      object.
#'      Any other values (e.g. \code{FALSE}) suppress the
#'      reference model calculations.
#'
#' @name log_results_to_excel
#' @return No return value.
#' @export
log_results_to_excel <- function(object,
                                 file,
                                 referencemodel = FALSE) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package \"openxlsx\" needed for this function to work.
          Please install it.",
         call. = FALSE)
  }

  if ( !grepl("\\.xlsx$", file, fixed = FALSE) ) file <- paste0(file, ".xlsx")
  .wrap_add_worksheet <- function(wb,
                                  x,
                                  sheetname,
                                  name = NA,
                                  colNames = FALSE,
                                  ...) {
    openxlsx::addWorksheet(wb, sheetname)
    srow <- 1
    if ( !is.na(name) ) {
      openxlsx::writeData(wb, sheetname, name)
      srow <- 2
    }
    openxlsx::writeData(wb, sheetname, x, startRow = srow, ...)
    openxlsx::setColWidths(wb,
                           sheetname,
                           cols = 1:(ncol(x) + !is.na(colNames)),
                           widths = "auto")
  }

  pds <- parseddata_summary(object)
  ydesc <- if ( inherits(pds[[1]], "character") ) {
    data.frame(summary_measure = names(pds[[1]]),
               outcome = pds[[1]])
  } else {
    as.data.frame(pds[[1]])
  }
  hps <- selected_hyperparameters(object)
  ces <- coefficients_summary(object)

  # production classification + hyperparameters:
  prod <- summary(object$prod$performance)[, -3]
  prod.hps <- data.frame(Measure = c("...",
                                      "Production Model Hyperparameters",
                                      "lambda",
                                      "alpha"),
                          Value = NA,
                          stringsAsFactors = FALSE)
  prod.hps[3, 2] <- object$prod$tuning$best$lambda
  prod.hps[4, 2] <- as.numeric(object$prod$tuning$best$alpha)
  prod <- as.data.frame(data.table::rbindlist(list(prod, prod.hps)))

  # do we need to calculate the reference models?:
  if ( identical(referencemodel, TRUE) ) {
    referencemodel <- refunreg(object)
  }

  if ( "refunreg" %in% class(referencemodel) ) {
    ces <- data.frame(ces, `...` = "-", coef_refunreg(referencemodel),
                      stringsAsFactors = FALSE)
  }

  # merge bits:
  sheets <- list(
    coversheet = setNames(as.data.frame(
      utils::capture.output(summary(object)),
      stringsAsFactors = FALSE
    ), ""),
    dqs.1 = ydesc,
    dqs.2 = as.data.frame(pds[[2]], stringsAsFactors = FALSE),
    classif = report_performance_summary(object),
    subclass = casesummary.performance(object$performance),
    hyperparameters = hps$joint,
    coefficients = ces,
    production = prod
  )

  cNames <- c(FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE)
  rNames <- c(FALSE, FALSE,  TRUE, FALSE, FALSE,  TRUE,  TRUE, FALSE)

  labs <- setNames(
    as.data.frame(
      rbind(
        c("Summary", NA),
        c("OutcomeDesc", "Descriptives for Outcome Variable"),
        c("PredictorDesc", "Descriptives for Matrix of Predictors"),
        c("CV-Performance",
          "Outer loop cross-validated performance"),
        c("SubjectClass", "Subject performance over outer loops"),
        c("Hyperparameters",
          "Hyperparameters chosen for outer loop folds/reps"),
        c("Coefficients", "Descriptives of Model coefficients"),
        c("Production Model", "Performance of Model (not cross-validated)")
      ), stringsAsFactors = FALSE),
    c("sheetname", "name")
  )

  cat(paste0("Opening: ", file, "\n"))
  wb <- openxlsx::createWorkbook()

  for (i in seq_along(sheets)) {
    .wrap_add_worksheet(wb,
                        x = sheets[[i]],
                        sheetname = labs$sheetname[i],
                        name = labs$name[i],
                        colNames = cNames[i],
                        rowNames = rNames[i])
  }

  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  cat(paste0("Write completed\n"))
  return(invisible(NULL))
}
