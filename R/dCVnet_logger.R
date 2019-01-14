
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
#'      reference models, or pass a pre-calculated \code{\link{reflogreg}}
#'      object.
#'      Any other values (e.g. \code{FALSE}) suppress the
#'      reference model caculations.
#'
#' @name log_results_to_excel
#' @return No return value.
#' @export
log_results_to_excel <- function(object,
                                 file,
                                 referencemodel = T) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package \"openxlsx\" needed for this function to work.
          Please install it.",
         call. = FALSE)
  }

  if ( !grepl("\\.xlsx$", file, fixed = F) ) { file <- paste0(file, ".xlsx") }
  .wrap_add_worksheet <- function(wb, x, sheetname, name = NA, colNames = F, ...) {
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
  hps <- selected_hyperparameters(object)
  ces <- coefficients_summary(object)

  # final classification + hyperparameters:
  final <- summary(classperformance(object$final))[,-3]
  #final <- data.frame(Measure = row.names(final), final)
  final.hps <- data.frame(Measure = c("...",
                                      "Final Model Hyperparameters",
                                      "lambda",
                                      "alpha"),
                          Value = NA,
                          stringsAsFactors = F)
  final.hps[3, 2] <- object$final$tuning$inner_best$lambda
  final.hps[4, 2] <- as.numeric(object$final$tuning$inner_best$alpha)
  final <- rbind(final, final.hps)

  # do we need to calculate the reference models?:
  if ( identical(referencemodel, TRUE) ) {
    referencemodel <- reflogreg(object)
  }

  if ( "reflogreg" %in% class(referencemodel) ) {
    ces <- data.frame(ces, `...` = "-", coef_reflogreg(referencemodel))
  }

  # merge bits:
  sheets <- list(
    coversheet = setNames(as.data.frame(
      utils::capture.output(summary(object))
    ),""),
    dqs.1 = setNames(as.data.frame(pds[[1]]),""),
    dqs.2 = as.data.frame(pds[[2]]),
    classif = report_classperformance_summary(object),
    subclass = casesummary.classperformance(object$performance),
    hyperparameters = hps$joint,
    coefficients = ces,
    final = final
  )

  cNames <- c(F, F, T, T, T, T, T, T)
  rNames <- c(F, F, T, F, F, T, T, F)

  labs <- setNames(
    as.data.frame(
      rbind(
        c("Summary", NA),
        c("OutcomeDesc", "Descriptives for Outcome Variable"),
        c("PredictorDesc", "Descriptives for Matrix of Predictors"),
        c("CV-Classification",
          "Outer loop cross-validated classification performance"),
        c("SubjectClass", "Subject classification over outer loops"),
        c("Hyperparameters", "Hyperparameters chosen for outer loop folds/reps"),
        c("Coefficients", "Descriptives of Model coefficients"),
        c("Final Model", "Performance of Final Model (not cross-validated)")
      ), stringsAsFactors = F),
    c("sheetname", "name")
  )

  cat(paste0("Opening: ", file, "\n"))
  wb <- openxlsx::createWorkbook()

  for (i in 1:length(sheets)) {
    .wrap_add_worksheet(wb,
                        x = sheets[[i]],
                        sheetname = labs$sheetname[i],
                        name = labs$name[i],
                        colNames = cNames[i],
                        rowNames = rNames[i])
  }

  openxlsx::saveWorkbook(wb, file = file, overwrite = T)
  cat(paste0("Write completed\n"))
  return(invisible(NULL))
}


