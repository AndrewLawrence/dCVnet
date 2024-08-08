# This script generates the blrsim dataset
#   used in the dCVnet-Worked-Example

set.seed(42)

# simulation settings:
nvar <- 40  # n predictors
ngood <- 10 # first ngood used for prediction
nobs_train <- 50   # observations in training dataset
nobs_test <- 1000  # observations in testing dataset
nobs_total <- nobs_train + nobs_test

# make a covariance matrix:
xsig <- diag(nvar)
# low background correlation:
xsig[upper.tri(xsig)] <- rnorm(n = (nvar * (nvar - 1)) / 2,
                               mean = 0.1, sd = 0.05)
# 3 highly correlated variable pairs:
xsig[1, 4] <- xsig[12, 13] <- xsig[13, 14] <- 0.8

# force symmetric positive definite:
xsig <- as.matrix(
  Matrix::nearPD(
    Matrix::forceSymmetric(xsig,
                           uplo = "U"),
    corr = TRUE
  )$mat
)

x <- MASS::mvrnorm(n = nobs_total,
                   mu = rep(0, nvar),
                   Sigma = xsig,
                   empirical = TRUE)
colnames(x) <- paste0("V", 1:nvar)

y <- ifelse(rowSums(x[, seq_len(ngood)]) > 0 +
              rnorm(n = nobs_total, mean = 0, sd = 2.5),
            yes = 1,
            no = 0)

# factoring out the error, how well would a classifier that knew the correct
#   rule do?
bayes.rate <- mean(y == as.numeric(rowSums(x[, 1:10]) > 0))
bayes.rate

# split off a training set of the first 25 of each outcome:
train_indices <- sort(c(which(y == 1)[1:25], which(y == 0)[1:25]))

x.test <- x[-train_indices, ]
y.test <- y[-train_indices]

x <- x[train_indices, ]
y <- y[train_indices]

blrsim <- list(x = x,
               y = y,
               x.test = x.test,
               y.test = y.test,
               bayes.rate = bayes.rate)

usethis::use_data(blrsim, overwrite = TRUE)

