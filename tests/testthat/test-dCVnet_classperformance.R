

context("tests of model performance objects")

# class performance object (from scratch)

# Setup
perfect_classification <- structure(
  data.frame(
    reference = c("A", "A", "B", "B"),
    prediction = c(0.15, 0.45, 0.65, 0.95),
    classification = c("A", "A", "B", "B"),
    label = c("example"),
    rowid = paste0("s", 1:4)
  ), class = c("performance", "data.frame")
)

imperfect_classification <- perfect_classification
imperfect_classification$reference <- rev(perfect_classification$reference)

perfect_classification.s <- summary(perfect_classification)
imperfect_classification.s <- summary(imperfect_classification)


perfect_classification.c <- suppressWarnings(
  calibration(perfect_classification))
imperfect_classification.c <- suppressWarnings(
  calibration(imperfect_classification))

class_measures <- c("Accuracy", "Kappa",
                    "Sensitivity", "Specificity",
                    "Pos Pred Value", "Neg Pred Value",
                    "Precision", "Recall",
                    "F1", "Prevalence", "Detection Rate",
                    "Detection Prevalence", "Balanced Accuracy",
                    "AUROC", "Brier")

test_that("classes are correct", {
  expect_s3_class(perfect_classification, "performance")
  expect_s3_class(imperfect_classification, "performance")
  expect_s3_class(perfect_classification.c, "calcoefs")
  expect_s3_class(imperfect_classification.c, "calcoefs")
})

# The tests:
test_that("binomial classification is as expected", {
  # perfect:
  expect_equal(perfect_classification.s$Value[
    perfect_classification.s$Measure %in% class_measures],
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 1, 1, 0.0875)
  )
  # imperfect:
  expect_equal(imperfect_classification.s$Value[
    imperfect_classification.s$Measure %in% class_measures],
    c(0, -1, 0, 0, 0, 0, 0, 0, NaN, 0.5, 0.0, 0.5, 0, 0, 0.5875)
  )
})

known_cal <- c(-11.62047, 55.67869)
test_that("binomial calibration is as expected", {
  # perfect:
  expect_equal(as.vector(perfect_classification.c), known_cal)
  # imperfect:
  expect_equal(as.vector(imperfect_classification.c), known_cal * -1)
})
