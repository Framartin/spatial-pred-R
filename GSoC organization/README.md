# GSoC organization

## 1st progress report

12/06/2015

For now, I have only worked on spdep.
- reading papers and doc
- code for in-sample predictors and one leave-one-out predictor of C.Thomas and al.
- add these ones to predict.sarlm()
- first update of the doc (still incomplete)
- some small optimizations
- (lot of) questions

## 2nd progress report

23/06/2015

- continue reading papers
- reorganize the code of predict.sarlm(): check, data, predictor. To avoid duplicated code. This allows the bug fix of NA in newdata to be apply for all models.
- add BP*, TS1, TC predictors
- add listw.decompose()
- repo reorganization + add README
- sarlm.examples.R (need to be updated)
- bug fixes

## Agenda 1st Skype call

26/06/2015
Solutions are under questions.

* class of the returned object: deal with in-sample / out-of-sample predictions ; deal with future prediction intervals (not now) / keep trend and signal attributes only for default predictor? 
  + If we choose to change the behavior, allow the use of the option legacy to use the old code. 
* How we ensure that the order of the spatial units in listw is the same than newdata?
  + see all methods for creating listw objects. 
* Is it possible to reorder spatial unit in a listw object? If not, I think we can work with only W in SparseMatrix, but it will imply modifications.
  + Easy solution: convert to Sparse Matrix and subset&order. But not efficient. For other solution: see spChk option and functions of spChkOption.R
* No Trend-Signal decomposition for the SAC model?
  + No
* What is aliased coef?
  + same definition than lm(): QR decomposition determine if X is full rank. If case of perfect collinearity, one variable is dropped, ie. marked as aliased. See predict.lm() to see what to do.
* Prevision case (ie. newdata for in-sample spatial units): should we use in-sample predictors? And "strange" cases where we have both in-sample and out-of-sample spatial units?
  + See what others functions do (glm, ts, etc). Use in-sample predictors. Print warning in the case of this "strange" case. But need to allow this case because it can happens often. In this case, use out-of-sample.
* Support log-transformation of the dependent variable?
  + see predict.lm(). For now, even if the log-transformation is quite popular in spatial econometrics, we keep it as a note for the future. I will implement predictors of Kato (2012) only if we have extra time at the end. Look at the total impact for prediction with X_s +1 to test if results are what is intended.

## 3rd progress report

03/07/2015

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

