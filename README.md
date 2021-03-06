# Project

This repository is a alpha version of the new `predict.sarlm()` function, which will be included into `spdep` in R-forge later.

The goal is to extend prediction in spatial econometrics for in-sample and out-of-sample spatial units by implementing predictors proposed by papers of the References section below. For more detail, see this [wiki page](https://github.com/rstats-gsoc/gsoc2015/wiki/Predict-methods-for-spatial-econometrics).

## Test

Run `library(spdep) ; source("sarlm.R") ; source("nb2mat.R")`. 

Examples are given in the `sarlm.examples.R` file.

## Documentation

The new documentation is provided in `predict.sarlm.Rd`.

## References

1. [Thomas-Agnan C, Laurent T, Goulard M (2014)](http://idei.fr/display.php?a=27788&site=TSE&data=TSE&lang=en). "About predictions in spatial autoregressive models: Optimal and almost optimal strategies", TSE Working Paper, n. 13-452, December 18, 2013, revised September 2014.
2. Bivand RS (2002). "Spatial Econometrics Functions in R: Classes and Methods." Journal of Geographical Systems, 4, 405-421.
3. [spatial-pred-r](http://htmlpreview.github.io/?https://github.com/jsay/spatial-pred-R/blob/master/DOC.html#sec-7) project
4. Kato, T. (2012). "Prediction in the lognormal regression model with spatial error dependence". Journal of Housing Economics 21: 66-76.
5.  Kelejian, H. H. and I. R. Prucha (2007). "The relative efficiencies of various predictors in spatial econometric models containing spatial lags". Regional Science and Urban Economics 37: 363-374.
6. Kato (2008)
7. ~~Griffith (2010) Spatial filtering and missing georeferenced data imputation: a comparison of the Getis and Griffith methods, in Perpectives on spatial data analysis, L. Anselin and S.J. Rey (eds), Springer-Verlag~~
8. LeSage and Pace (2004)


Notes on the previous references:

1. Definition of in-sample and out-of-sample prediction in the spatial case. Using its notations (eg. for the sub-spatial weight matrix). Proposed new predictors for the LAG model. Extended here to SDM model (and to SAC model, with warning ; to be confirmed).
2. Description of the previous version of `predict.sarlm()`. The `pred.sarlm` class may have to change.
3. Using its notations for new predictors. Using a different approach for predictors: for `mixed` models we include WX in the X matrix object. Our goal is to include this custom function `sppred()` into the `predict.sarlm()` framework.
4. Prediction with a log-transformed variable on SEM model. Anti-log transformation is biased: proposal of 4 other predictors to correct (partially or completely) this bias
5. 3 optimal predictors for the SAC model, based on 3 different information sets + 2 intuitive predictors (with one for the error model)
6. TODO: BPn predictor for SEM model?
7. ~~EM procedure for missing value with Spatial Filtering~~
8. TODO: BP predictor in-sample for SEM model


## Other useful References

* Bivand RS, Piras G (2015) "Comparing Implementations of Estimation Methods for Spatial Econometrics." Journal of Statistical Software, 63(18), 1-36
* Millo G, Piras G (2012). "splm: Spatial Panel Data Models in R." Journal of Statistical Software, 47(1), 1-38.
* Millo G (2014). "Maximum likelihood estimation of spatially and serially correlated panels with random effects". Computational Statistics & Data Analysis 71 (March), 914-933.
* Piras G (2010). "sphet: Spatial Models with Heteroskedastic Innovations in R." Journal of Statistical Software, 35(1), 1-21.


## License

The following files are distributed under the GPL (>= 2) license ([GPLv2](https://cran.r-project.org/web/licenses/GPL-2) and [GPLv3](https://cran.r-project.org/web/licenses/GPL-3)):

- `sarlm.R`
- `sarlm.examples.R`
- `nb2mat.R`
- `predict.sarlm.Rd`
- `nb2mat.Rd`


## Notes

This project is partially based on the work of Jean-Sauveur AY, Raja CHAKIR and Julie LE GALLO
