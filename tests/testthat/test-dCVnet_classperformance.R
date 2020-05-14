
# class performance object (from scratch)

context("tests of model performance objects")

perfect_classification <- structure(
  data.frame(
    reference = c("A", "A", "B", "B"),
    prediction = c(0.15, 0.45, 0.65, 0.95),
    classification = c("A", "A", "B", "B"),
    label = c("example"),
    rowid = paste0("s", 1:4)
  ), class = c("performance", "data.frame")
)

perfect_classification.s <- summary(perfect_classification)

class_measures <- c("Accuracy", "Kappa",
                    "Sensitivity", "Specificity",
                    "Pos Pred Value", "Neg Pred Value",
                    "Precision", "Recall",
                    "F1", "Prevalence", "Detection Rate",
                    "Detection Prevalence", "Balanced Accuracy", "AUROC")


test_that("prediction measures work", {
  expect_equal(perfect_classification.s$Value[
    perfect_classification.s$Measure %in% class_measures],
               c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 1, 1)
  )
})
