# GSoC organization

## 1st progress report

For now, I have only worked on spdep.
- reading papers and doc
- code for in-sample predictors and one leave-one-out predictor of C.Thomas and al.
- add these ones to predict.sarlm()
- first update of the doc (still incomplete)
- some small optimizations
- (lot of) questions

## 2nd progress report

- continue reading papers
- reorganize the code of predict.sarlm(): check, data, predictor. To avoid duplicated code. This allows the bug fix of NA in newdata to be apply for all models.
- add BP*, TS1, TC predictors
- add listw.decompose()
- repo reorganization + add README
- sarlm.examples.R (need to be updated)
- bug fixes

## Agenda 1st Skype call

* class of the returned object: deal with in-sample / out-of-sample predictions ; deal with future prediction intervals (not now) / keep trend and signal attributes only for default predictor?
* How we ensure that the order of the spatial unit in listw is the same than newdata?
* Is it possible to reorder spatial unit in a listw object? If not, I think we can work with only W in SparseMatrix, but it will imply modifications.
* No Trend-Signal decomposition for the SAC model?
* What is aliased coef?
* Prevision case (ie. newdata for in-sample spatial units): should we use in-sample predictors? And "strange" cases where we have both in-sample and out-of-sample spatial units?
* Support log-transformation of the dependent variable?

## 3rd progress report

* Update sarlm.examples.R
* Fix some bugs
* Add BP predictor
* Confirmation from C.Thomas that predictors of the LAG model can be extend to SDM model by replacing X_S by (X_S, WX_S), X_O by (X_O, WX_O) and beta by (beta, teta)' where teta is the vector of coefficients associated with WX_S (already done in sarlm.R)

## TODO

* contact C. Thomas: questions sent. Waiting for her advices on SEM model, code for EM approach, and corrections of her paper
* merge code from spatial-pred-r (need adaptations)
* reorder listw (or SparseMatrix?)
* WARNING: newdata != prediction of out-of-sample spatial units. How to deal with the difference between out-of-sample prediction (in the sense of C.Thomas & al. (2015)) and prevision: both have newdata. Does the second case need in-sample predictors?
* warnings if spatial units are on newdata and data
* add spatial units names for rows
* work on splm & spdep
