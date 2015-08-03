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
  + See how other functions handle that. If we choose to change the behavior, allow the use of the option legacy to use the old code. 
* How we ensure that the order of the spatial units in listw is the same than newdata?
  + see all methods for creating listw objects. 
* Is it possible to reorder spatial unit in a listw object? If not, I think we can work with only W in SparseMatrix, but it will imply modifications.
  + Easy solution: convert to Sparse Matrix and subset&order. But not efficient. For other solution: see spChk option and functions of spChkOption.R (e.g chkIDs())
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

## 4th progress report

23/07/2015

* Read ref 7
* Important clean of the code of sarlm.R (including listw issues previously described, re-organization, clean checks, fix TS not defined for out-of-sample, etc.)
* Add TC1/KP1 & trend predictor types
* Add KP4 alias for TS1 predictor type (same definition)
* Add spChk option
* Update examples
* Add error if newdata as no rownames
* Add call attribute
* Add support of matrix/Matrix argument for invIrW()
* Fix minor bugs

## 5th progress report

* Checks and fix issues with W and listw for out-of-sample
* Add spChk option for out-of-sample
* Clean some old TODOs
* add region.id attribute in sarlm.pred class which is passed to rownames in as.data.frame.sarlm.pred()
* warning if spatial units are on newdata and data
* add the prevision case for some in-sample predictor

## TODO

* contact C. Thomas: questions sent. Waiting for her advices on SEM model, code for EM approach, and corrections of her paper
* merge code from spatial-pred-r (need adaptations)
* WARNING: newdata != prediction of out-of-sample spatial units. How to deal with the difference between out-of-sample prediction (in the sense of C.Thomas et al. (2015)) and prevision: both have newdata. Does the second case need in-sample predictors?
* work on splm & spdep
* update invIrW() documentation
* IMPORTANT BUG: Error of computation of WXo: shloud we compute WoXo or out-of-sample row of WX? Should we keep the ex behavior, ie WoXo? Only with legacy=TRUE?
* Read chapter on BLUP of the Encyclopedia of GIS, Springer (2007)
