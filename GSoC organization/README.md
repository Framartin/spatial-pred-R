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

## TODO

* contact C. Thomas
* merge code from spatial-pred-r (need adaptations)
* reorder listw (or SparseMatrix?)
