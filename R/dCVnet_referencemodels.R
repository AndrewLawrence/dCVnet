
# Reference Models --------------------------------------------------------


reference_models <- function(x, ...) {
  UseMethod("reference_models", x)
}


reference_models.dCVnet <- function(object) {
  parsed <- object$input$parsed
  n <- min(table(parsed$y))
  # above we assume effective N for a logistic regression is n minority cases.
  k <- ncol(parsed$x_mat)
  ideal_nvars <- round(n/5) # ideally we want at least 5 cases per predictor.

  if ( k > ideal_nvars ) {
    # if we have more variables than subjects then first
    #   reduce the number of variables via a (svd)PCA on the data.
    x_pca <- prcomp(parsed$x_mat,
                    rank. = ideal_nvars,
                    retx = T)
    pglm.data <- data.frame(y = parsed$y, as.data.frame.matrix(x_pca$x))

    pglm <- glm(y ~ . , data = pglm.data, family = "binomial")
  } else {
    pglm.data <- data.frame(y = parsed$y, parsed$x_mat)
    pglm <- glm(y ~ ., data = pglm.data, family = "binomial")
  }
  # next univariate prediction:
  punivariate <- lapply(1:k, function(i) {
    f <- as.formula(paste0("y ~ ", colnames(parsed$x_mat)[i]))
    glm(f, data = data.frame(y = parsed$y, parsed$x_mat), family = "binomial")
  } )
  names(punivariate) <- colnames(parsed$x_mat)
  return(list(glm = pglm,
              univariate = structure(punivariate,
                                     class = c("glmlist",
                                               class(punivariate)))))
}


report_reference_classperformance_summary <- function(refobj) {

  glm <- summary(classperformance(refobj$glm), "GLM")

  glm <- glm[,-3]
  names(glm)[2] <- "Reference GLM"

  univ <- report_classperformance_summary(blarg1.ref$univariate)

  names(univ)[2:5] <- paste("UnivPred", names(univ)[2:5])

  return(data.frame(glm, univ[,-1]))
}

