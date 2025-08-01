# ~ -----------------------------------------------------------------------
# Data info:
#   https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.info.txt
# ~ -----------------------------------------------------------------------

# Prostate data info
#
# Predictors (columns 1--8)
#
#   lcavol
#   lweight
#   age
#   lbph
#   svi
#   lcp
#   gleason
#   pgg45
#
# outcome (column 9)
#
#   lpsa
#
# train/test indicator (column 10)
#
# This last column indicates which 67 observations were used as the
# "training set" and which 30 as the test set, as described on page 48
# in the book.
#
# There was an error in these data in the first edition of this
# book. Subject 32 had a value of 6.1 for lweight, which translates to a
# 449 gm prostate! The correct value is 44.9 gm. We are grateful to
# Prof. Stephen W. Link for alerting us to this error.
#
# The features must first be scaled to have mean zero and  variance 96 (=n)
# before the analyses in Tables 3.1 and beyond.
# That is, if x is the  96 by 8 matrix of features, then
#   we compute: xp <- scale(x,TRUE,TRUE)


# ~ -----------------------------------------------------------------------

# import the data
# (original source:
#     https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data)
prostate <- read.table("data-raw/prostate.data")

# write it out:
usethis::use_data(prostate, overwrite = TRUE)

# ~ -----------------------------------------------------------------------
