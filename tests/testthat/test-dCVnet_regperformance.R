# regression-type performance objects (from scratch) --------------------------

# make some fake data for testing:
data <- list(
  perfect_preds = structure(
    data.frame(
      reference = 1:10,
      prediction = 1:10,
      label = "perfect"
    ),
    class = c("performance", "data.frame"),
    family = "gaussian"
  ),
  mean_preds = structure(
    data.frame(
      reference = 1:10,
      prediction = rep(5.5, 10),
      label = "mean-only"
    ),
    class = c("performance", "data.frame"),
    family = "gaussian"
  ),
  offset_mean_preds = structure(
    data.frame(
      reference = 1:10,
      prediction = rep(0, 10),
      label = "offset-mean"
    ),
    class = c("performance", "data.frame"),
    family = "gaussian"
  ),
  offset_perfect_preds = structure(
    data.frame(
      reference = 1:10,
      prediction = (1:10) - 3,
      label = "perfect"
    ),
    class = c("performance", "data.frame"),
    family = "gaussian"
  ),
  degenerate_preds = structure(
    data.frame(
      reference = 0,
      prediction = 2,
      label = "degenerate"
    ),
    class = c("performance", "data.frame"),
    family = "gaussian"
  )
)


# known-good values:
known_good <- list(
  perfect_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE"
      ),
      Value = c(0, 0, 1, 0, 1, 0, 0, 0),
      label = c(
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -8L)
  ),
  mean_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE"
      ),
      Value = c(
        2.87228132326901,
        2.5,
        0,
        0,
        NA,
        8.25,
        0.948683298050514,
        0.82572282384477
      ),
      label = c(
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -8L)
  ),
  offset_mean_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE"
      ),
      Value = c(
        6.20483682299543,
        5.5,
        -3.66666666666667,
        5.5,
        NA,
        38.5,
        2.04939015319192,
        1.81659021245849
      ),
      label = c(
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -8L)
  ),
  offset_perfect_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE"
      ),
      Value = c(
        3,
        3,
        -0.0909090909090908,
        3,
        1,
        9,
        0.990867388613724,
        0.990867388613724
      ),
      label = c(
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -8L)
  ),
  degenerate_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE"
      ),
      Value = c(2, 2, NaN, NA, NA, 4, NA, NA),
      label = c(
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -8L)
  )
)

test_that("no errors during calculations", {
  expect_no_error({
    summaries <<- lapply(data, summary)
  })
})


test_that("Calibration function rejects regression input:", {
  expect_error(calibration(data[[1]]))
})

test_that("calculated values are OK:",
          expect_equal(summaries, known_good))



# as above, but poisson ---------------------------------------------------

pdata <- list(
  perfect_preds = structure(
    data.frame(
      reference = 1:10,
      prediction = 1:10,
      label = "perfect"
    ),
    class = c("performance", "data.frame"),
    family = "poisson"
  ),
  mean_preds = structure(
    data.frame(
      reference = 1:10,
      prediction = rep(5.5, 10),
      label = "mean-only"
    ),
    class = c("performance", "data.frame"),
    family = "poisson"
  ),
  offset_mean_preds = structure(
    data.frame(
      reference = 1:10,
      prediction = rep(1, 10),
      label = "offset-mean"
    ),
    class = c("performance", "data.frame"),
    family = "poisson"
  ),
  offset_perfect_preds = structure(
    data.frame(
      reference = 1:10,
      prediction = (1:10) + 3,
      label = "perfect"
    ),
    class = c("performance", "data.frame"),
    family = "poisson"
  ),
  degenerate_preds = structure(
    data.frame(
      reference = 0,
      prediction = 2,
      label = "degenerate"
    ),
    class = c("performance", "data.frame"),
    family = "poisson"
  )
)


pknown_good <- list(
  perfect_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE",
        "poisson_r2",
        "RMSLE",
        "poisson_cal_Intercept",
        "poisson_cal_Slope"
      ),
      Value = c(0, 0, 1, 0, 1, 0, 0, 0, 1, 0,
                -0.666312196652238, 1.24996422652761),
      label = c(
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -12L)
  ),
  mean_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE",
        "poisson_r2",
        "RMSLE",
        "poisson_cal_Intercept",
        "poisson_cal_Slope"
      ),
      Value = c(
        2.87228132326901,
        2.5,
        0,
        0,
        NA,
        8.25,
        0.948683298050514,
        0.82572282384477,
        0,
        0.539894766372882,
        -0.167054083383958,
        NA
      ),
      label = c(
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only",
        "mean-only"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -12L)
  ),
  offset_mean_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE",
        "poisson_r2",
        "RMSLE",
        "poisson_cal_Intercept",
        "poisson_cal_Slope"
      ),
      Value = c(
        5.33853912601566,
        4.5,
        -2.45454545454545,
        4.5,
        NA,
        28.5,
        1.76326144093578,
        1.48630108292059,
        -5.85952751966043,
        1.18073388247787,
        1.01160091295769,
        NA
      ),
      label = c(
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean",
        "offset-mean"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -12L)
  ),
  offset_perfect_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE",
        "poisson_r2",
        "RMSLE",
        "poisson_cal_Intercept",
        "poisson_cal_Slope"
      ),
      Value = c(
        3,
        3,
        -0.0909090909090908,
        -3,
        1,
        9,
        0.990867388613724,
        0.990867388613724,
        0.107960165946634,
        0.49567144652133,
        -2.64718789520931,
        1.89962638329521
      ),
      label = c(
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect",
        "perfect"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -12L)
  ),
  degenerate_preds = structure(
    list(
      Measure = c(
        "RMSE",
        "MAE",
        "r2",
        "cal_Intercept",
        "cal_Slope",
        "Brier",
        "SDScaledRMSE",
        "SDScaledMAE",
        "poisson_r2",
        "RMSLE",
        "poisson_cal_Intercept",
        "poisson_cal_Slope"
      ),
      Value = c(2, 2, NaN, NA, NA, 4, NA, NA, -Inf, 1.09861228866811, NA, NA),
      label = c(
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate",
        "degenerate"
      )
    ),
    class = "data.frame",
    row.names = c(NA, -12L)
  )
)

test_that("no errors during poisson calculations", {
  expect_no_error({
    psummaries <<- lapply(pdata, summary)
  })
})

test_that("Calibration function rejects poisson-reg input:",
          {
            expect_error(calibration(pdata[[1]]))
          })

test_that("calculated values are OK:",
          expect_equal(psummaries, pknown_good))
