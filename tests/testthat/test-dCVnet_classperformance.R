
# class performance object (from scratch)

context("tests of model performance objects")

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

class_measures <- c("Accuracy", "Kappa",
                    "Sensitivity", "Specificity",
                    "Pos Pred Value", "Neg Pred Value",
                    "Precision", "Recall",
                    "F1", "Prevalence", "Detection Rate",
                    "Detection Prevalence", "Balanced Accuracy",
                    "AUROC", "Brier")

# The tests:
test_that("perfect binomial classification is as expected", {
  expect_equal(perfect_classification.s$Value[
    perfect_classification.s$Measure %in% class_measures],
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 1, 1, 0.0875)
  )
})

test_that("imperfect binomial classification is as expected", {
  expect_equal(imperfect_classification.s$Value[
    imperfect_classification.s$Measure %in% class_measures],
    c(0, -1, 0, 0, 0, 0, 0, 0, NaN, 0.5, 0.0, 0.5, 0, 0, 0.5875)
  )
})
