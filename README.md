# Dev version

This is the dev version of dCVnet. Currently being updated to support
other model families than binomial. Documentation may be out of date and
functionality only partly supported.

# dCVnet

**dCVnet** is an R package to estimate **d**oubly
**C**ross-**V**alidated Elastic-**net** regularised generalised linear
models (glm) with an approachable interface. dCVnet adds nested repeated
k-fold cross-validation, and convenience functions, to the regularised
glm fit by the [glmnet
package](https://cran.r-project.org/web/packages/glmnet/index.html).

If you use dCVnet in your research please cite our paper [Lawrence et al
(2021)](https://pubmed.ncbi.nlm.nih.gov/34175478/)

## Getting Started

A working installation of [R](https://www.r-project.org/) is required. A
fully featured interface to R (e.g. [RStudio](https://www.rstudio.com/))
is recommended.

dCVnet is not (yet) on [CRAN](https://cran.r-project.org/), so the
remotes package is useful to download and build from github:

#### First, install remotes (if needed):

    install.packages("remotes")

The commands below will install missing package dependencies (see the
[DESCRIPTION](DESCRIPTION) file *Imports* section). It will then run a
toy example from the package's main function.

#### Install dCVnet (from GitHub):

    remotes::install_github("AndrewLawrence/dCVnet", dependencies = TRUE, build_vignettes = TRUE)

#### -or- install the dev version of dCVnet (from GitHub):

    remotes::install_github("AndrewLawrence/dCVnet@dev", dependencies = TRUE, build_vignettes = TRUE)

#### Install dCVnet (from an Archive):

    remotes::install_local("path/to/dCVnet_1.0.8.tar.gz", dependencies = TRUE, build_vignettes = TRUE)

#### Run a simple example:

    library(dCVnet)
    example(dCVnet, run.dontrun = TRUE)
    # to see the usage without running calculations set run.dontrun = FALSE

#### List of dCVnet Vignettes

    vignette(package = "dCVnet")

Note: this needs `build_vignettes = TRUE` to be set at installation.

## Support

Please search for your issue in the project
[Issues](https://github.com/AndrewLawrence/dCVnet/issues) section. If
that doesn't clear it up please make a new issue, if possible try to
give a Reproducible Example (see
[here](https://stackoverflow.com/a/5963610), or
[here](http://adv-r.had.co.nz/Reproducibility.html)).

## Notes

-   This package is not affiliated with glmnet or it's authors.
-   [AndrewLawrence](https://github.com/AndrewLawrence)'s work on this
    software is funded by the National Institute for Health Research
    (NIHR) Biomedical Research Centre at South London and Maudsley NHS
    Foundation Trust and King's College London. The views expressed are
    those of the author(s) and not necessarily those of the NIHR or the
    Department of Health and Social Care.

## License

This project is licensed under the
[GPL\>3](https://www.gnu.org/licenses/gpl.html). See
[DESCRIPTION](DESCRIPTION) file.

## Authors

-   **Andrew Lawrence** -
    [AndrewLawrence](https://github.com/AndrewLawrence)

# What is dCVnet?

[ see the Presentations folder for slides from a talk on dCVnet given
2021-09-15 ]

The motivating problem behind dCVnet is prediction
modelling<sup>[1](#fn1)</sup> in data with relatively few observations
(*n*) for the number of predictors (*p*), especially where there may be
uninformative or redundant predictors which the analyst isn't willing,
or able, to remove.

In an ideal world we would collect more observations (i.e. increase
*n*), or better understand which predictors to include/exclude or how to
best model them (i.e. reduce *p*), but this can be impractical or
impossible. For example, it is often a necessary step to justify funding
for further research to increase *n* or reduce *p*.

With few observations and many predictors several inter-related
statistical problems arise. These problems become
worse<sup>[2](#fn2)</sup> with greater ratios of *p*/*n*:

-   [Overfitting](https://en.wikipedia.org/wiki/Overfitting)
-   [Multicollinearity](https://en.wikipedia.org/wiki/Multicollinearity)
-   [need for Variable
    Selection](https://en.wikipedia.org/wiki/Feature_selection)

*dCVnet* uses *elastic-net regularisation* (from glmnet) to combat these
problems. *double cross-validation*<sup>[3](#fn3)</sup> is applied to
tune the regularisation and validly assess model performance.

#### Cross-validation for generalisation error

A model which is overfit is tuned to the noise in the sample rather than
reproducible relationships. As a result it will perform poorly in new
(unseen) data. This failure to perform well in new data is termed
generalisation (or out-of-sample) error.

Generalisation error can be assessed using properly conducted
cross-validation. The model is repeatedly refit in different subsets of
the data and performance evaluated in the observations which were not
used for model fitting. Appropriately cross-validated estimates of model
performance are unaffected by the optimism caused by overfitting and
reflect the likely performance of the model in unseen data.

There are different forms of cross-validation, *dCVnet* implements
repeated k-fold cross-validation.

However, cross-validation only tells the analyst if overfitting is
occurring, it is not a means to reduce overfitting. For this purpose we
can apply **regularisation** which produces more cautious models which
are likely to generalise better to new data.

#### Elastic-net Regularisation

Regularisation adds a cost to the complexity of the model. Unregularised
models optimise the fit of the model to the data, regularised models
optimise the fit of the model given a budget of allowable complexity.
This results in shrinkage of model coefficients towards zero. This makes
models more cautious and can substantially improve generalisation to
unseen data. *dCVnet* uses elastic-net regularisation from the
[glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)
package for R.

The type of regularisation used by *dCVnet* ([Elastic-net
Regularisation](https://en.wikipedia.org/wiki/Elastic_net_regularization))
is a combination of two types of regularisation with the aim of avoiding
their weaknesses and benefiting from their strengths:

-   Ridge regression (using a *L2-*penalty) allows predictors with
    perfect *collinearity*, but every predictor contributes (the
    solution is not sparse).

-   LASSO (Least Absolute Shrinkage and Selection Operator) regression
    uses the *L1-*penalty. It produces variable selection effect by
    favouring a sparse solution (meaning less important coefficients
    drop to zero), however LASSO is unstable when working with
    correlated predictors.

#### Elastic-net Hyperparameters

Adding regularisation to a model introduces "algorithm
[hyperparameters](https://en.wikipedia.org/wiki/Hyperparameter_(machine_learning))" -
these are settings which which must be tuned/optimised for each problem.

Elastic-net regularisation requires two hyperparameters be specified:

-   `alpha` - the balance of *L1*- and *L2*-regularisation penalties.
    (L2/Ridge only : alpha = 0, L1/LASSO only : alpha = 1)
-   `lambda` - penalty factor determining the combined amount of
    regularisation.

There are no *default* values for these parameters, suitable values vary
depending on the problem and so should be 'tuned'.

One way to tune parameters without overfitting is to use
Cross-validation to select values which perform well in unseen data.
This is a form of model selection.

If the cross-validation for hyperparameter selection is combined with
the cross-validation to estimate generalisation error this will add back
in optimism to our estimates of the generalisation error.

To combat this cross-validation for generalisation error must be
completely independent of the cross-validation for hyperparameter
selection, see [Crawley & Talbot (2010; JMLR
11:2079-2107)](http://www.jmlr.org/papers/v11/cawley10a.html) for a
fuller description of the issue. Nesting the hyperparameter tuning can
achieve this.

#### Nested Cross-validation

Double cross-validation<sup>[4](#fn4)</sup> is implemented to allow
principled (and independent) selection of the optimal hyperparameters
for generalisability, and estimation of performance in out-of-sample
data when hyperparameters are estimated in this way. Double
cross-validation is computationally expensive, but ensures
hyperparameter tuning is fully separate from performance evaluation.

### Summary

This package aims to provide an approachable interface for conducting a
nested (or double) cross-validation of the elastic net solution to a
two-class prediction problem. The Elastic-net calculations (and some
inner loop cross-validation) are performed by the R package
[glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)
which dCVnet depends on.

dCVnet was developed to work with two-class outcomes (binary logistic),
but work on implementing other model families supported by glmnet is
ongoing.

### Footnotes

<a name="fn1"><sup>1</sup></a> dCVnet can be useful
for inference, but this is not its main purpose. The
time-consuming outer cross-validation loop is not as important for
inference, other software can be used directly.

<a name="fn2"><sup>2</sup></a> Where *p*/*n* \> 1, the standard
least-squares regression solutions are not defined and generalised
models will have convergence problems. In both cases predictors will
have perfect multicollinearity. Where *p*/*n* \> 0.1, common [rules of
thumb](https://en.wikipedia.org/wiki/One_in_ten_rule) for sample size
are violated.

<a name="fn3"><sup>3</sup></a> Double cross-validation is also called
*nested* or *nested-loop* cross-validation.

With less flexible models, and enough data the optimism which nested-CV
addresses can be negligible. However, nested cross-validation is
particularly important with smaller datasets. Demonstrating internal
validity without validation leakage is important for [reproducible
research](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5238707/).

<a name="fn4"><sup>4</sup></a> Other examples of nested CV in R:
[MLR](https://pat-s.github.io/mlr/articles/tutorial/devel/nested_resampling.html),
[TANDEM](https://www.rdocumentation.org/packages/TANDEM/versions/1.0.2/topics/nested.cv),
[nlcv](https://cran.r-project.org/web/packages/nlcv/vignettes/nlcv.pdf),
[caret/rsample](http://appliedpredictivemodeling.com/blog/2017/9/2/njdc83d01pzysvvlgik02t5qnaljnd)
[nestedcv](https://github.com/myles-lewis/nestedcv)
