
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