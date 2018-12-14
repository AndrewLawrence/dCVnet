
# TODO: Implement lambda + 1se, lambda + 3pc, lambda + Xse, lambda + Xpc.

# TODO: S3 class + print / summary methods for reference model objects.
#         - include clear description of whether PCA reduction was applied.
#
# TODO: Add 'tidy' methods for objects / summaries / results.
#
# TODO: Write something which dumps model summary and 'complete' results to
#         an excel file.
#
# TODO: Summary methods take a looong time to run when you have even a moderate
#         number of reps (~100). Can we optimise the code?
#
# TODO: Sampling options - Leave-group-out,
#                          Proportional sampling by outcome or covariate/
#       This could be achieved by switching to rsample?
#       UPDATE: 14/12/2018 - the current code is already proportional by outcome.
#                             could add option to stratify by other covariates.
#                            rsample seems to be overkill as it stores the data.
#
# TODO: add checks that randomised folds are unique?
#
# TODO: Timing values are incorrect when the model is fit in parallel. Nice to have a parallel progress bar.







# Fixed (05/12/2018):
# TODO: If you request only a single rep then the summary measures 'break' with errors or display issues.

# FIXED (29/11/2018):
# TODO: when running in parallel the inner-loop (per-alpha) report lines are confusing. Slowfix: investigate a parallel progress bar. Quick-fix: suppress these when parallel.
#   added an mc.cores check to the inner-loop print command.

# TODO: ROC averaging.
#   average_rocdata included.

# NOT REPLICATED (19/11/2018)
# TODO: On installation: "Error e1071 required" probably needs adding to the description/imports/depends.
# RESULTS: systematic search found no functions from e1071.
# Possibly the error message relates to another package.
# UPDATE: it was the confusionMatrix function from caret.
#         Caret only 'suggests' e1071, so it is not installed as a dependency.
#         adding e1071 to the imports fixes the error.


# FIXED (26/10/2018)
# Needs print functions for new S3 types.
# Needs better summary function for dCVnet

# FIXED (25/10/2018)
# Currently dCVnet:::predict_cat.glm() fails if
# the outcome is is not the first variable in the data.frame.

# FIXED (25/10/2018)
# Currently code fails if factor labels don't produce legal column names.

# OLD:
#
#
# Version 0.7.0
#
#   - Add balanced sampling options for train/test folds (inner and/or outer).
#       either  for y (class imbalance)
#               or for x (e.g. representative proportion of M/F)
#
#   - Principled averaging of ROC curves
#
#   - develop a plot/diagnostic which shows (per-subject) the outer-loop fold
#       performance for when that subject is test vs. when they are train.
#       This may help identify hard to predict participants.
#
#   - tests for / evaluation of data where p>>n.
#
#   - add report of performance within the inner CV,
#       (this may help with diagnosis of:
#           underfitting / overfitting / too few iterations)
#
#   - Standardisation of Variables.
#       currently the preprocessing (scaling) is performed at the level of the
#       outer loop before entering the inner loop.
#       Thus there is some theoretical leakage between the test/train data
#       for the inner loop.
#       This may detriment hyperparameter selection.
#       However, any fix will result in a additional processing overhead as
#       quadratically more calls to caret::preProcess will be required.
#       Investigate this and potentially change the implmentation.
#       UPDATE: 27/07/2018
#         standardise argument is not being set.
#         So inner loop is being standardised prior to running.
#         (returning coefficients in original scale.)
#         There is no leakage here, but standardisation happens 2x!
#
#   - make a some S3 classes:
#         dCVnet (print, predict, summary)
#         hyperparam.dCVnet (plot)
#         ROC.dCVnet (plot)
#     and convert existing log, plot functions to make use of methods.
#
#   - convert to package (in progress).
#
#   # Not needed:
#   - write up and submit bug to glmnet package authors.
#         It turns out this has been done:
# https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R