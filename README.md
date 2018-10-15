# dCVnet

This is the source code of the R Package 'dCVnet' which estimates
**d**oubly **C**ross-**V**alidated Elastic-**net** regularised logistic models.

The aim is to provide an approachable interface for conducting a nested (or double) 
cross-validation of the elastic net solution to a two-class prediction problem. 
Double cross-validation is computationally expensive, but ensures 
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
(see the [DESCRIPTION](DESCRIPTION) file *Imports* section). It will run a 
fast toy example from the package's main function.

In most cases nested cross-validation will be run on larger, more complex
problems than the given example. It will also likely need many more iterations
to give reliabile results.


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


