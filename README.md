# dCVnet

This is the source code of the R Package 'dCVnet' which estimates
**d**oubly **C**ross-**V**alidated Elastic-**net** regularised logistic models.

# Motivation

The motivating problem behind dCVnet is valid prediction modelling in datasets 
with relatively few observations (`n`) for the number of predictors (`p`). 
In typical modelling problems the predictors may be highly correlated and
include a mix of predictive and non-predictive variables.
Over-fitting is common or inevitable in these circumstances (e.g. if `p > n`). 
Elastic-net regularisation permits predictors with collinearity
and performs variable selection (through L1-penalisation - 
Least Absolute Shrinkage and Selection Operator: LASSO).

However two hyperparameters must be specified in elastic-net regularisation:
the balance of L1- and L2-penalties (`alpha`), and the total ammount of 
regularisation (`lambda`).

The double cross-validation approach allows principled and independent 
estimation of both the performance in out-of-sample data, and 
the optimal hyperparameters for generalisability.

This package aims to provide an approachable interface for conducting a nested 
(or double) cross-validation of the elastic net solution to a two-class 
prediction problem. Double cross-validation is computationally expensive, but ensures 
hyperparameter tuning is fully separate from performance evaluation. 
See [Crawley & Talbot JMLR 11(2010):2079-2107](http://www.jmlr.org/papers/v11/cawley10a.html) 
for a description of the issue.

The Elastic-net calculations (and some inner loop cross-validation) are 
performed by the R package
[glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) which this
package depends on.

The development of this package was motivated by a lack of options for 
double/nested cross-validation. Although it currently only works with two-class
outcomes a future development goal is to make it work with the other outcome 
families catered for by glmnet (gaussian, poisson, multinomial, cox, mgaussian).

## Getting Started

A working installation of [R](https://www.r-project.org/) is required.
This package is not yet on [CRAN](https://cran.r-project.org/), 
so currently the devtools package is helpful to easily download and build 
from source.

The commands below will install missing package dependencies 
(see the [DESCRIPTION](DESCRIPTION) file *Imports* section). It will then
run a toy example from the package's main function.

Note: for most real analyses problems will be larger and more iterations will be
required to give reliable results. dCVnet is computationally expensive.

Regularised regression methods are not helpful if the model is 


```
if ( ! require("devtools") ) install.packages("devtools")
devtools::install_github("AndrewLawrence/dCVnet", dependencies = TRUE)

library(dCVnet)

example(dCVnet, run.dontrun = T)
```

## Notes
* This package is very much in development and comes with 
ABSOLUTELY NO WARRANTY - please use at your own risk.


## License

This project is licensed under the 
[GPL>3](https://www.gnu.org/licenses/gpl.html) as set out in the 
[DESCRIPTION](DESCRIPTION) file.


## Authors

* **Andrew Lawrence** - [AndrewLawrence](https://github.com/AndrewLawrence)


