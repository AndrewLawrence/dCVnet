---
output:
  html_document: default
  pdf_document: default
---
# dCVnet

This is the source code of the R Package 'dCVnet' which estimates
**d**oubly **C**ross-**V**alidated Elastic-**net** regularised logistic models.

### Overview & Rationale

The motivating problem behind dCVnet is prediction modelling in datasets 
with relatively few observations (*n*)<sup>[1](#fn1)</sup> for the number of 
predictors (*p*), especially where there may be uninformative or redundant 
predictors which the analyst isn't willing, or able, to remove.

These circumstances lead to a number of related statistical problems which 
become worse with greater ratios of *p*/*n*:<sup>[2](#fn2)</sup>

* [Overfitting](https://en.wikipedia.org/wiki/Overfitting)
* [Multicollinearity](https://en.wikipedia.org/wiki/Multicollinearity)
* [Variable Selection](https://en.wikipedia.org/wiki/Feature_selection)

*dCVnet* implements *double cross-validation*<sup>[3](#fn3)</sup> and 
*elastic-net regularisation* to appropriately analyse data in these
circumstances.

Properly conducted 
[cross-validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)) 
gives estimates of model performance which are unaffected by the optimism 
inherent in overfitting. Cross-validated estimates better reflect performance 
for unseen data. *dCVnet* implements k-fold cross-validation (which is 
optionally repeated).
However, cross-validation itself can only tell you if overfitting is occurring, 
it is not a treatment for overfitting. More cautious models which 
generalise better to new data are required.

[Regularisation](https://en.wikipedia.org/wiki/Regularization_\(mathematics\)) 
adds a cost to the complexity of the model. This results in [shrinkage](https://en.wikipedia.org/wiki/Shrinkage_estimator) of 
model coefficients. This makes models more cautious and can substantially 
improve generalisation to unseen data. 
*dCVnet* uses elastic-net regularisation provided by the [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) package for
R.

Adding regularisation to a model introduces algorithm 
[hyperparameters](https://en.wikipedia.org/wiki/Hyperparameter_(machine_learning))
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
and performs *variable selection* (through *L1*-penalisation - 
the Least Absolute Shrinkage and Selection Operator; LASSO).

Elastic-net regularisation requires two hyperparameters be specified:

* `alpha` - the balance of *L1*- and *L2*-regularisation penalties.
* `lambda` - the total amount of regularisation.

there are no default values for these parameters, and suitable values 
vary depending on the problem.

Double cross-validation<sup>[4](#fn4)</sup> is implemented in *dCVnet*
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
double/nested cross-validation. Although it currently only works with two-class
outcomes a future development goal is to implement for other outcome 
families catered for by glmnet (`gaussian`, `poisson`, `multinomial`, `cox`,
`mgaussian`).

## Getting Started

A working installation of [R](https://www.r-project.org/) is required.
This package is not on [CRAN](https://cran.r-project.org/), 
so the devtools package is useful to download and build 
from source.

The commands below will install missing package dependencies 
(see the [DESCRIPTION](DESCRIPTION) file *Imports* section). It will then
run a toy example from the package's main function.


```
if ( ! require("devtools") ) install.packages("devtools")
devtools::install_github("AndrewLawrence/dCVnet", dependencies = TRUE)

library(dCVnet)

example(dCVnet, run.dontrun = T)
```

## Notes
* This package is not affiliated with the glmnet package.
* This package is under active development and comes with 
ABSOLUTELY NO WARRANTY - use at your own risk.

## License

This project is licensed under the 
[GPL>3](https://www.gnu.org/licenses/gpl.html) as set out in the 
[DESCRIPTION](DESCRIPTION) file.


## Authors

* **Andrew Lawrence** - [AndrewLawrence](https://github.com/AndrewLawrence)


<br><br>
<hr>

<a name="fn1">1</a>: For two-class prediction problems 
(e.g. binary logistic regression) the class-balance impacts sample size, 
the rarer of the two outcomes should be used.


<a name="fn2">2</a>: Where *p*/*n* > 1 standard least-squares solutions 
(for continuous outcomes) are not defined. For binary logistic models 
convergence problems are likely. In both cases predictors will have
perfect mutlicollinearity. Where *p*/*n* > 0.1,
conventional [rules of thumb](https://en.wikipedia.org/wiki/One_in_ten_rule) 
for logistic regression sample size are violated.


<a name="fn3">3</a>: Double cross-validation is alternatively known as *nested* 
or *nested-loop* cross-validation.

<a name="fn4">4</a>: Other examples of nested CV in R:
[MLR](https://pat-s.github.io/mlr/articles/tutorial/devel/nested_resampling.html),
[TANDEM](https://www.rdocumentation.org/packages/TANDEM/versions/1.0.2/topics/nested.cv),
[nlcv](https://cran.r-project.org/web/packages/nlcv/vignettes/nlcv.pdf),
[caret/rsample](http://appliedpredictivemodeling.com/blog/2017/9/2/njdc83d01pzysvvlgik02t5qnaljnd)
