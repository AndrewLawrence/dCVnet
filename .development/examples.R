library(glmnet)
library(dCVnet)

#      )                              (
#   ( /(                    )         )\ )
#   )\())  (  (     (    ( /(    (   (()/(  (      ) (  (
#  ((_)\  ))\ )(   ))\   )\())  ))\   /(_)) )(  ( /( )\))( (   (    (
#   _((_)/((_|()\ /((_) ((_)\  /((_) (_))_ (()\ )(_)|(_))\ )\  )\ ) )\
#  | || (_))  ((_|_))   | |(_)(_))    |   \ ((_|(_)_ (()(_|(_)_(_/(((_)
#  | __ / -_)| '_/ -_)  | '_ \/ -_)   | |) | '_/ _` / _` / _ \ ' \)|_-<
#  |_||_\___||_| \___|  |_.__/\___|   |___/|_| \__,_\__, \___/_||_|/__/
#                                                   |___/


# Example Application -----------------------------------------------------

# Use iris data:
data(BinomialExample, package = "glmnet") # x and y

df <- data.frame(y = factor(y,
                            levels = c(0,1),
                            labels = c("control", "case")),
                 x)

# Parse data into dCVnet input format:
parsed <- dCVnet::parse_dCVnet_input(f = y ~ .,
                                     data = df,
                                     positive = "case")

# a simple cvglmnet:
m1 <- cv.glmnet(x = parsed$x_mat, y = parsed$y, type.measure = "deviance",
                family = "binomial", alpha = 0.5)
plot(m1)
coef(m1)
table(m1 = factor(predict(m1,
                          newx = parsed$x_mat,
                          s = "lambda.min",
                          type = "class"),
                  levels = levels(parsed$y)),
      actual = parsed$y)

m2 <- cv.glmnet(x = parsed$x_mat, y = parsed$y, type.measure = "deviance",
                family = "binomial", alpha = 0.5, grouped = F)


# ~ Single alpha examples ---------------------------------------------------

# set up fixed folds:
fixfolds <- lapply(1:7, function(i) {
  caret::createFolds(y = parsed$y, k = 10, list = FALSE, returnTrain = FALSE)
})

alphalist <- c(0.2, 0.5, 0.8)


maxlambdas <- lambda_rangefinder(x = parsed$x_mat,
                                 y = parsed$y,
                                 alphalist = alphalist,
                                 prop = 0.9)

lambdalist <- lambda_seq_list(maxlambdas, 1000, 0.0001)


test1 <- dCVnet:::repeated.cv.glmnet(fixfolds,
                                    x = parsed$x_mat,
                                    y = parsed$y,
                                    type.measure = "class",
                                    lambdas = lambdalist[[2]],
                                    alpha = 0.5)

test2 <- dCVnet:::repeated.cv.glmnet(fixfolds,
                                    x = parsed$x_mat,
                                    y = parsed$y,
                                    type.measure = "deviance",
                                    lambdas = lambdalist[[2]],
                                    alpha = 0.5)

test3 <- dCVnet:::repeated.cv.glmnet(fixfolds,
                                    x = parsed$x_mat,
                                    y = parsed$y,
                                    type.measure = "auc",
                                    lambdas = lambdalist[[2]],
                                    alpha = 0.5)



plot(test1) # (mis)class(ification)
summary(test1)

plot(test2) # deviance
summary(test2)

plot(test3) # auc
summary(test3)

# ~ Multialpha examples: ---------------------------------------------------

# for multialpha we now have a
tmp1 <- dCVnet:::multialpha.repeated.cv.glmnet(alphalist = c(0.2, 0.5, 0.8),
                                               lambdas = lambdalist,
                                               k = 5, nrep = 5,
                                               x = parsed$x_mat,
                                               y = parsed$y,
                                               type.measure = "class")

tmp2 <- dCVnet:::multialpha.repeated.cv.glmnet(alphalist = c(0.2, 0.5, 0.8),
                                               lambdas = lambdalist,
                                               k = 5, nrep = 5,
                                               x = parsed$x_mat,
                                               y = parsed$y,
                                               type.measure = "deviance")

tmp3 <- dCVnet:::multialpha.repeated.cv.glmnet(alphalist = c(0.2, 0.5, 0.8),
                                               lambdas = lambdalist,
                                               k = 5, nrep = 5,
                                               x = parsed$x_mat,
                                               y = parsed$y,
                                               type.measure = "auc")

plot(tmp1)
summary(tmp1)

plot(tmp2)
summary(tmp2)

plot(tmp3)
summary(tmp3)


# ~ Outer loop example ------------------------------------------------------

# We now also need to specify outer loop cv stuff as well as
#   empirical cutoff and the parameters for
#     lambda_rangefinder and lambda_seq_list
blarg1 <- dCVnet::dCVnet(f = y ~ .,
                         data = df,
                         positive = "case",
                         alphalist = c(0.2, 0.5, 0.8),
                         k_inner = 10, nrep_inner = 3,
                         k_outer = 5, nrep_outer = 2,
                         option.empirical_cutoff = FALSE,
                         type.measure = "class",
                         nlambda = 100)

blarg2 <- dCVnet::dCVnet(f = y ~ .,
                         data = df,
                         positive = "case",
                         alphalist = c(0.2, 0.5, 0.8),
                         k_inner = 10, nrep_inner = 3,
                         k_outer = 10, nrep_outer = 2,
                         option.empirical_cutoff = FALSE,
                         type.measure = "auc",
                         nlambda = 100)

str(blarg1,1)

# we can also plot all the inner results together:
plot(blarg1)
plot(blarg2)

# We can still get access to individual inner tuning plots:
plot(blarg1$tuning[[1]]$tuning) # rep1
plot(blarg1$tuning[[2]]$tuning) # rep2
plot(blarg1$tuning[[3]]$tuning) # rep3

# We can extract and visualise the per-fold coefficients:
coef(blarg1, type = "mean")
plot_outerloop_coefs(blarg1)
plot_outerloop_coefs(blarg2)

# And classification stats / AUROC:
summary(classperformance(blarg1))
summary(classperformance(blarg2))

blarg1.ref <- dCVnet_refmodels(blarg1)

summary(classperformance(blarg1.ref$glm), "GLM")
report_classperformance_summary(blarg1.ref$univariate)

report_reference_classperformance_summary(blarg1.ref)

# any classperformance can be turned into a rocplot.
head(extract_rocdata(classperformance(blarg1.ref$univariate)))

plot(extract_rocdata(classperformance(blarg1)), legend = T)
plot(extract_rocdata(classperformance(blarg1$final)))
plot(extract_rocdata(classperformance(blarg1.ref$glm)))
plot(extract_rocdata(classperformance(blarg1.ref$univariate)), legend = T)

# HERE.
# TODO: nice formatting for numbers / writing outputs to excel / pdf.
# TODO: parallel processing, n.cores?
# TODO: average of ROC curves?
# TODO: per-subject model performances?

# Example per-subject:
blarg1$performance %>%
  group_by(rowid) %>%
  summarise(mean = mean(classification == reference)) %>%
  select(mean) %>%
  table() %>%
  as.data.frame() %>% setNames(c("ProportionCorrect", "NSubjects")) %>%
  ggplot(aes(y = NSubjects, x = ProportionCorrect)) +
  geom_bar(stat = "identity", fill = "goldenrod", colour = "black") +
  geom_hline(yintercept = 0.0) +
  xlab("Proportion of reps subject correctly classified") +
  theme_light()

# What should the object summary tell us...

summary(tmp1)
summary(tmp1, print = F)

summary(blarg1)


summary(blarg1$tuning$OutFold1.Rep1$tuning)
