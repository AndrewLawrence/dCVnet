# This script generates the modsel_perf results
#   used in the dCVnet-limitations vignette
library(dCVnet)
set.seed(42)
y <- rep(0:1, length.out = 150)
# 10,000 predictors:
X <- data.frame(matrix(rnorm(150 * 10000), nrow = 150, ncol = 10000))

# Screening for 20 best predictors by R2:
tests <- apply(X, 2,
               function(x) {
                 summary(lm(y ~ ., data = data.frame(y = y, x = x)))$r.squared
               })
Xbest <- sort(tests, decreasing = TRUE)[1:20]
# X6096 X3530 X7954  X180  X527 X6679  X259 X3112  X739 X7395 X8452 X6370
# 0.111 0.085 0.078 0.078 0.075 0.071 0.069 0.066 0.066 0.065 0.062 0.061
# X7124 X7741 X4644 X6175 X4216 X9230 X7945 X5058
# 0.061 0.060 0.060 0.060 0.060 0.059 0.059 0.059

# Cheating:
#   (Use the 20 'best' predictors)
dn1 <- dCVnet(y = y, data = X[, names(X) %in% names(Xbest)], nrep_outer = 10)

# Playing Fair:
#   (Use all predictors)
dn2 <- dCVnet(y = y, data = X, nrep_outer = 10)

# extract performance:
modsel_perf <- list(
  Cheating = report_performance_summary(dn1),
  PlayingFair = report_performance_summary(dn2)
)

# report mean, minimum and maximum Accuracy and AUROC:
modsel_perf <- lapply(modsel_perf, function(x) {
  subset(x,
         subset = x$Measure %in% c("Accuracy", "AUROC"),
         select = c("Measure", "mean", "min", "max"))
})

usethis::use_data(modsel_perf, overwrite = TRUE)
