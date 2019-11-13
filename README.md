---
output:
  html_document: default
  pdf_document:
    highlight: tango
urlcolor: blue
---
# dCVnet

**dCVnet** is an R package to estimate 
**d**oubly **C**ross-**V**alidated Elastic-**net** 
regularised generalised linear models (glm). dCVnet adds nested 
cross-validation, and convenience functions, to regularised glm fit by the 
[glmnet package](https://cran.r-project.org/web/packages/glmnet/index.html).

#### Overview & Rationale

The motivating problem behind dCVnet is prediction modelling^[1](#fn1)^ in data 
with relatively few observations (*n*) for the number of 
predictors (*p*), especially where there may be uninformative or redundant 
predictors which the analyst isn't willing, or able, to remove.

In an ideal world we would collect more observations (i.e. increase *n*), 
or better understand which predictors to include/exclude or how to best model 
them (i.e. reduce *p*), but this can be impractical/expensive.
It is often a necessary step to justify funding for further research to
increase *n* or reduce *p*. An analysis of albeit limited data can help 
to support this.

With few observations and many predictors several inter-related 
statistical problems arise and these become worse^[2](#fn2)^ with greater ratios
of *p*/*n*:

* [Overfitting](https://en.wikipedia.org/wiki/Overfitting)
* [Multicollinearity](https://en.wikipedia.org/wiki/Multicollinearity)
* [Variable Selection](https://en.wikipedia.org/wiki/Feature_selection)

*dCVnet* implements *elastic-net regularisation* to combat these problems 
with *double cross-validation*^[3](#fn3)^ to tune the regularisation and 
validly assess model performance.

#### Cross-validation for generalisation error

A model which is overfit is tuned to the noise in the sample rather than 
reproducible relationships. As a result it will perform poorly in new (unseen) 
data. This failure to perform well in new data is termed 
[generalisation (or out-of-sample) error](https://en.wikipedia.org/wiki/Generalization_error).

Generalisation error can be assessed using properly conducted internal- or 
[cross-validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)). 
In cross-validation the model is repeatedly refit in different subsets of the 
data and performance evaluated in the cases which were not used for model 
fitting. Cross-validated estimates of model performance are unaffected by the 
optimism caused by overfitting and reflect the likely performance of the model 
in unseen data.

There are different forms of cross-validation, *dCVnet* implements repeated 
k-fold cross-validation.

However, cross-validation only tells the analyst if overfitting is occurring, 
it cannot reduce overfitting. For this purpose we use **regularisation** which 
produces more cautious models likely to generalise better to new data.

#### Elastic-net Regularisation

[Regularisation](https://en.wikipedia.org/wiki/Regularization_\(mathematics\)) 
adds a cost to the complexity of the model. Unregularised models optimise the 
fit of the model to the data, regularised models optimise the fit of the model
given a budget of allowable complexity. This results in [shrinkage](https://en.wikipedia.org/wiki/Shrinkage_estimator) of 
model coefficients towards zero. This makes models more cautious and can 
substantially improve generalisation to unseen data. 
*dCVnet* uses elastic-net regularisation from the  [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) package for
R.

The type of regularisation employed by *dCVnet* ([Elastic-net](https://en.wikipedia.org/wiki/Elastic_net_regularization))
allows predictors with perfect *collinearity*
and performs *variable selection* (by incorporating *L1*-penalisation - 
also called the Least Absolute Shrinkage and Selection Operator; LASSO).

#### Elastic-net Hyperparameters

However, adding regularisation to a model introduces "algorithm 
[hyperparameters](https://en.wikipedia.org/wiki/Hyperparameter_(machine_learning))"
- these are settings which which must be tuned/optimised for each problem.

Elastic-net regularisation requires two hyperparameters be specified:

* `alpha` - the balance of *L1*- and *L2*-regularisation penalties.
(L2 only : alpha = 0, L1 only : alpha = 1)
* `lambda` - the total amount of regularisation.

there are no *default* values for these parameters, and suitable values 
vary depending on the problem. Some considerations are:

* Problems with a greater number of relevant variables (and/or stronger relationships for those relevant predictors) will be fit best with a greater lambda.
* Problems which are sparse (i.e. many predictors are irrelevant) are better fit with more L1-like models.

One way to tune parameters without overfitting is to use Cross-validation 
to select values which perform well in unseen data. This is a form of model 
selection.

If the cross-validation for hyperparameter selection is combined with the 
cross-validation to estimate generalisation error this will add back in 
optimism to our estimates of the generalisation error. 

To combat this cross-validation for generalisation error must be completely 
independent of the cross-validation for hyperparameter selection, see 
[Crawley & Talbot (2010; JMLR 11:2079-2107)](http://www.jmlr.org/papers/v11/cawley10a.html) 
for a fuller description of the issue. Nesting the hyperparameter tuning 
can achieve this.


#### Nested Cross-validation

Double cross-validation^[4](#fn4)^ is implemented to allow principled (and independent) selection 
of the optimal hyperparameters for generalisability, and estimation of 
performance in out-of-sample data when hyperparameters are estimated in this 
way. Double cross-validation is computationally expensive, 
but ensures hyperparameter tuning is fully separate from performance evaluation. 

### Summary

This package aims to provide an approachable interface for conducting a nested 
(or double) cross-validation of the elastic net solution to a two-class 
prediction problem. The Elastic-net calculations (and some inner loop cross-validation) are 
performed by the R package
[glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) which this
package depends on.

The development of this package was motivated by a lack of options for 
double/nested cross-validation. dCVnet currently works with 
the binary logistic model for two-class outcomes, but in the future other 
outcome families in glmnet 
(`gaussian`, `poisson`, `multinomial`, `cox`, `mgaussian`) will be implemented.

## Getting Started

A working installation of [R](https://www.r-project.org/) is required.
A fully featured interface to R (e.g. [RStudio](https://www.rstudio.com/)) 
is recommended.

dCVnet is not (yet) on [CRAN](https://cran.r-project.org/),
so the devtools package is useful to download and build 
from source:

#### First, install devtools:
```
install.packages("devtools")
```

The commands below will install missing package dependencies 
(see the [DESCRIPTION](DESCRIPTION) file *Imports* section). It will then
run a toy example from the package's main function.

#### Install dCVnet from Github:
```
devtools::install_github("AndrewLawrence/dCVnet", dependencies = TRUE)
```

#### Install dCVnet from Archive:
```
devtools::install_local("path/to/dCVnet_1.0.2.tar.gz", dependencies = TRUE)
```

#### Run a simple example:
```
library(dCVnet)
example(dCVnet, run.dontrun = TRUE)
# to see the usage without running calculations set run.dontrun = FALSE
```

## Notes
* This package is not affiliated with glmnet or it's authors.
* This package is under active development and comes with ABSOLUTELY NO WARRANTY 
- use at your own risk.
* This software represents independent research part funded by the National Institute for Health Research (NIHR) Biomedical Research Centre at South London and Maudsley NHS Foundation Trust and King’s College London. The views expressed are those of the author(s) and not necessarily those of the NIHR or the Department of Health and Social Care.

## License

This project is licensed under the 
[GPL>3](https://www.gnu.org/licenses/gpl.html) as set out in the 
[DESCRIPTION](DESCRIPTION) file.


## Authors

* **Andrew Lawrence** - [AndrewLawrence](https://github.com/AndrewLawrence)


\newpage

### Footnotes

<a name="fn1">^1^</a> regularised-models/dCVnet can be useful for inference, 
but this is not really their main purpose, and the time-consuming outer 
cross-validation loop is not very important for inference.

<a name="fn2">^2^</a> Where *p*/*n* > 1, the standard least-squares regression 
solutions are not defined and generalised models will have convergence problems.
In both cases predictors will have perfect mutlicollinearity. Where *p*/*n* > 0.1,
common [rules of thumb](https://en.wikipedia.org/wiki/One_in_ten_rule) 
for sample size are violated.

<a name="fn3">^3^</a> Double cross-validation is also called *nested* 
or *nested-loop* cross-validation.

With stable models with enough data the optimism which nested-CV addresses can 
be negligible. This has lead some to conclude that nested-cv is 
not very important to do nested-CV (e.g. Max Kuhn, author of 
[caret](http://topepo.github.io/caret/index.html), in this
[blog post](http://appliedpredictivemodeling.com/blog/2017/9/2/njdc83d01pzysvvlgik02t5qnaljnd)). 
However, nested cross-validation is particularly important with smaller datasets 
and unstable models and demonstrating internal validity without validation leakage is
important for reproducible research. [Publication guidelines](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5238707/) are moving to address this.


<a name="fn4">^4^</a> Other examples of nested CV in R:
[MLR](https://pat-s.github.io/mlr/articles/tutorial/devel/nested_resampling.html),
[TANDEM](https://www.rdocumentation.org/packages/TANDEM/versions/1.0.2/topics/nested.cv),
[nlcv](https://cran.r-project.org/web/packages/nlcv/vignettes/nlcv.pdf),
[caret/rsample](http://appliedpredictivemodeling.com/blog/2017/9/2/njdc83d01pzysvvlgik02t5qnaljnd)
