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

04/08/2015

* Checks and fix issues with W and listw for out-of-sample
* Add spChk option for out-of-sample
* Clean some old TODOs
* add region.id attribute in sarlm.pred class which is passed to rownames in as.data.frame.sarlm.pred()
* warning if spatial units are on newdata and data
* minor changes

## 6th progress report

14/08/2015

* fix computation of lagged variables (introduce new way of computing WXo by default and legacy.mixed option to keep the old behavior)
* add the forecast case for some in-sample predictor (newdata has the same rownames than object$y)
* add KP2 and KP3 predictors
* add KP5 for SEM model (for a SAC, warning)
* documentation update
* remove warnings for TC, TS and TS1 predictors when used on a sac model
* TS1 with a error model now returns trend (bug fix)
* minor changes

## 7th progress report

24/08/2015

* add BPW, BPN predictors
* fix typo
* update documentation (predict.sarlm.Rd)
* update examples
* clean repo

## Final fixes

* Add license (GPL >= 2)
* merged modifs from CRAN by R.Bivand
* add BP1, BPW1 and BPN1 types
* update and improve documentation
* important bug fix in xx1 predictors: Wi need to be re-normalized in order to return the same results if we execute predict.sarlm() independantly or with multiple rows on newdata


## TODO

* correct bug if newdata has only one row (with drop=F when subseting elements). Done?
* work on splm ~~& sphet~~ (not enough scientific paper for sphet)
* update predict.sarlm() documentation (include a summary table)
* support all.data option for other predictors
* add XWy from spatial-pred-r? rename it?
* For leave-one-out predictors, if an out-of-sample unit has no neighbors, warns that non-spatial predictor will do the same more efficiently (?)
* support listw2 for sac model (problem: not developed by Kelejian et al (2007))
* add pred.se for KPx predictors
* Read chapter on BLUP of the Encyclopedia of GIS, Springer (2007)

