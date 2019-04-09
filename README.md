---
output:
  html_document: default
  pdf_document:
    highlight: tango
urlcolor: blue
---
# dCVnet

This is the source code of the R Package 'dCVnet' which estimates
**d**oubly **C**ross-**V**alidated Elastic-**net** regularised logistic models.

### Overview & Rationale

The motivating problem behind dCVnet is prediction modelling in datasets 
with relatively few observations (*n*)^[1](#fn1)^ for the number of 
predictors (*p*), especially where there may be uninformative or redundant 
predictors which the analyst isn't willing, or able, to remove.

In an ideal world we would collect more data, but this is often impractical
or too expensive.

The described circumstances produce several interelated problems which 
become worse^[2](#fn2)^ with greater ratios of *p*/*n*.

* [Overfitting](https://en.wikipedia.org/wiki/Overfitting)
* [Multicollinearity](https://en.wikipedia.org/wiki/Multicollinearity)
* [Variable Selection](https://en.wikipedia.org/wiki/Feature_selection)

*dCVnet* implements *double cross-validation*^[3](#fn3)^ and 
*elastic-net regularisation* to combat these problems.

Properly conducted 
[cross-validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)) 
gives estimates of model performance unaffected by the optimism 
caused by overfitting. Cross-validated estimates reflect the likely performance 
of the model in unseen data.
*dCVnet* implements repeated k-fold cross-validation.
However, cross-validation only tells the analyst if overfitting is occurring, 
it cannot correct overfitting. Regularisation produces more cautious models 
which will generalise better to new data.

[Regularisation](https://en.wikipedia.org/wiki/Regularization_\(mathematics\)) 
adds a cost to the complexity of the model. This results in [shrinkage](https://en.wikipedia.org/wiki/Shrinkage_estimator) of 
model coefficients towards zero. This makes models more cautious and can 
substantially improve generalisation to unseen data. 
*dCVnet* uses elastic-net regularisation from the  [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) package for
R.

Adding regularisation to a model introduces "algorithm 
[hyperparameters](https://en.wikipedia.org/wiki/Hyperparameter_(machine_learning))"
- settings which which must be tuned/optimised to each problem.
Cross-validation allows this tuning to select hyperparameters which produce 
models that give the best performance for unseen data.
However, to avoid optimism bias, cross-validation for hyperparameter selection
must be completely independent of the cross-validation which estimates model 
performance. 
See 
[Crawley & Talbot JMLR 11(2010):2079-2107](http://www.jmlr.org/papers/v11/cawley10a.html) 
for a description of the issue.

The type of regularisation employed by *dCVnet* ([Elastic-net](https://en.wikipedia.org/wiki/Elastic_net_regularization))
permits predictors with perfect *collinearity*
and performs *variable selection* (by incorporating *L1*-penalisation - 
also called the Least Absolute Shrinkage and Selection Operator; LASSO).

Elastic-net regularisation requires two hyperparameters be specified:

* `alpha` - the balance of *L1*- and *L2*-regularisation penalties.
* `lambda` - the total amount of regularisation.

there are no default values for these parameters, and suitable values 
vary depending on the problem.

Double cross-validation^[4](#fn4)^ is implemented in *dCVnet*
to allow principled and independent selection 
of the optimal hyperparameters for generalisability, and estimation of 
performance in out-of-sample data when hyperparameters are estimated in this 
way. Double cross-validation is computationally expensive, 
but ensures hyperparameter tuning is fully separate from performance evaluation. 

This package aims to provide an approachable interface for conducting a nested 
(or double) cross-validation of the elastic net solution to a two-class 
prediction problem. The Elastic-net calculations (and some inner loop cross-validation) are 
performed by the R package
[glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) which this
package depends on.

The development of this package was motivated by a lack of options for 
double/nested cross-validation. Although it currently only works with 
the binary logistic model for two-class outcomes the aim is to eventually 
implement the other outcome families in glmnet 
(`gaussian`, `poisson`, `multinomial`, `cox`, `mgaussian`).

## Getting Started

A working installation of [R](https://www.r-project.org/) is required.
This package is not on [CRAN](https://cran.r-project.org/), 
so the devtools package is useful to download and build 
from source.

The commands below will install missing package dependencies 
(see the [DESCRIPTION](DESCRIPTION) file *Imports* section). It will then
run a toy example from the package's main function.

#### from Github
```
if ( ! require("devtools") ) install.packages("devtools")
devtools::install_github("AndrewLawrence/dCVnet", dependencies = TRUE)
```

#### from Archive
```
if ( ! require("devtools") ) install.packages("devtools")
devtools::install_local("path/to/dCVnet_1.0.2.tar.gz", dependencies = TRUE)
```

#### Quickstart example:
```
library(dCVnet)

example(dCVnet, run.dontrun = T)
```

## Notes
* This package is not affiliated with the glmnet package.
* This package is under active development and comes with 
ABSOLUTELY NO WARRANTY - use at your own risk.
* This software represents independent research part funded by the National Institute for Health Research (NIHR) Biomedical Research Centre at South London and Maudsley NHS Foundation Trust and King’s College London. The views expressed are those of the author(s) and not necessarily those of the NIHR or the Department of Health and Social Care.

## License

This project is licensed under the 
[GPL>3](https://www.gnu.org/licenses/gpl.html) as set out in the 
[DESCRIPTION](DESCRIPTION) file.


## Authors

* **Andrew Lawrence** - [AndrewLawrence](https://github.com/AndrewLawrence)


\newpage

<a name="fn1">^1^</a> For two-class prediction problems like binary logistic
regression any class-imbalance has an impact on the effective sample size. 
*n* should be number of the rarer of the two outcomes.


<a name="fn2">^2^</a> Where *p*/*n* > 1, the standard least-squares solutions 
for continuous outcomes are not defined. For generalised models 
convergence problems are likely. In both cases predictors will have
perfect mutlicollinearity. Where *p*/*n* > 0.1,
common [rules of thumb](https://en.wikipedia.org/wiki/One_in_ten_rule) 
for logistic regression sample size are violated.


<a name="fn3">^3^</a> Double cross-validation is also called *nested* 
or *nested-loop* cross-validation.

<a name="fn4">^4^</a> Other examples of nested CV in R:
[MLR](https://pat-s.github.io/mlr/articles/tutorial/devel/nested_resampling.html),
[TANDEM](https://www.rdocumentation.org/packages/TANDEM/versions/1.0.2/topics/nested.cv),
[nlcv](https://cran.r-project.org/web/packages/nlcv/vignettes/nlcv.pdf),
[caret/rsample](http://appliedpredictivemodeling.com/blog/2017/9/2/njdc83d01pzysvvlgik02t5qnaljnd)
