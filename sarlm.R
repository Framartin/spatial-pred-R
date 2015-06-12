# Copyright 2002-12 by Roger Bivand
#

residuals.sarlm <- function(object, ...) {
	if (is.null(object$na.action))
		object$residuals
	else napredict(object$na.action, object$residuals)
}

deviance.sarlm <- function(object, ...) {
	object$SSE
}

coef.sarlm <- function(object, ...) {
	ret <- NULL
#	ret <- sqrt(object$s2)
#	names(ret) <- "sigma"
	if (object$type == "error") ret <- c(ret, object$lambda)
	else if (object$type == "lag" || object$type == "mixed")
            ret <- c(ret, object$rho)
        else if (object$type == "sac" || object$type == "sacmixed")
            ret <- c(ret, object$rho, object$lambda)
	ret <- c(ret, object$coefficients)

	ret
}

vcov.sarlm <- function(object, ...) {
	if (object$ase) res <- object$resvar[-1,-1]
        else {
            if (!is.null(object$fdHess)) {
                if (object$insert) res <- object$resvar[-1,-1]
                else res <- object$resvar
            } else {
                stop("vcov not available for this model")
            }
        }
        res
}


fitted.sarlm <- function(object, ...) {
	if (is.null(object$na.action))
		object$fitted.values
	else napredict(object$na.action, object$fitted.values)
}


# retourne la valeur de la prédiction + un attributs trend et signal
predict.sarlm <- function(object, newdata=NULL, listw=NULL, type=NULL, 
                          zero.policy=NULL, legacy=TRUE, power=NULL, order=250, tol=.Machine$double.eps^(3/5), #pred.se=FALSE, lagImpact=NULL, 
                          ...) {
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (object$type == "sac") stop("no predict method for sac") # Wy et We
  if (is.null(power)) power <- object$method != "eigen"
  stopifnot(is.logical(legacy))
  stopifnot(is.logical(power))
  #        if (pred.se && object$type == "error") {
  #            pred.se <- FALSE
  #            warning("standard error estimates not available for error models")
  #        }
  #        if (pred.se && is.null(lagImpact))
  #            stop("lagImpact object from impact method required for standard error estimate")
  if (is.null(newdata)) { # in-sample pred
    X <- object$X
    B <- object$coefficients
    y <- object$y
    tarX <- object$tarX
    tary <- object$tary
    trend <- X %*% B
    if (is.null(type)) { # defaut predictor
      res <- fitted.values(object)
      if (object$type == "error") { # We ou WX+We
        attr(res, "trend") <- as.vector(trend)
        attr(res, "signal") <- as.vector( -1 * (tary - y) - 					-1 * (tarX - X) %*% B)
      } else { # lag model or other
        attr(res, "trend") <- as.vector(trend)
        attr(res, "signal") <- as.vector( -1 * (tary - y))
      }
    } else { # new predictors
      if (! type %in% c("X", "TC", "BP")) stop("no such predictor type")
      #TODO: lag model only? + adapt to SDM model
      if (is.null(listw) || !inherits(listw, "listw"))
        stop ("spatial weights list required")
      if (nrow(X) != length(listw$neighbours))
        stop("mismatch between data and spatial weights")
      if (type == "X") { # warning: equals to WX predictor for mixed model
        res <- as.vector(trend)
      } else if (type %in% c("TC", "BP")) { # need to compute TC
        if (power){
          W <- as(listw, "CsparseMatrix")
          TC <- powerWeights(W, rho = object$rho, X = X, order = order, tol = tol) %*% B # TODO: not good for SDM model!
        } else {
          TC <- invIrW(listw, object$rho) %*% trend # TODO: not good for SDM model!
        }
      }
      if(type == "TC") res <- as.vector(TC)
      if(type == "BP") {
        W <- as(listw, "CsparseMatrix")
        Qss <- 1/object$s2 * (Diagonal(dim(W)[1]) - (object$rho * t(W))) %*% (Diagonal(dim(W)[1]) - (object$rho * W)) # precision matrix for LAG model
        DiagQss <- Diagonal(x = diag(Qss))
        BP <- TC - solve(DiagQss) %*% (Qss - DiagQss) %*% (y - TC)
        # Can BP also be applied to the SEM model? Cf LeSage and Pace (2004). Note: \hat{\mu_i} need to be adapted
        res <- as.vector(BP)
      }
      attr(res, "trend") <- as.vector(trend)
      attr(res, "signal") <- NULL
    }
  }
  else { # out-of-sample
    if (is.null(type)) {
      if (object$type == "error") {
        if (object$etype == "error") { # We
          B <- object$coefficients
          #			tt <- terms(object$lm.model) 
          #			X <- model.matrix(delete.response(tt), data=newdata)
          frm <- formula(object$call)
          mt <- delete.response(terms(frm, data=newdata)) # returns a terms object for the same model but with no response variable
          #			mf <- lm(object$formula, newdata, method="model.frame")
          mf <- model.frame(mt, newdata)
          X <- model.matrix(mt, mf)
          #  accommodate aliased coefficients 120314     # TODO : WHAT ?
          if (any(object$aliased)) # TODO : What ?
            X <- X[,-which(object$aliased)]
          trend <- X %*% B
          signal <- rep(0, length(trend))
          res <- trend + signal # pourquoi ajouter signal si c'est un vecteur de 0 ?
          attr(res, "trend") <- trend
          attr(res, "signal") <- signal
        } else if (object$etype == "emixed") { # WX + We
          if (is.null(listw) || !inherits(listw, "listw")) 
            stop ("spatial weights list required")
          if (nrow(newdata) != length(listw$neighbours)) # TODO: need to be changed for new predictors
            stop("mismatch between newdata and spatial weights")
          B <- object$coefficients
          #			mt <- terms(object$formula, data = newdata)
          frm <- formula(object$call)
          mt <- delete.response(terms(frm, data=newdata))
          #			mf <- lm(object$formula, newdata, method="model.frame")
          mf <- model.frame(mt, newdata)
          X <- model.matrix(mt, mf)
          K <- ifelse(colnames(X)[1] == "(Intercept)", 2, 1)
          m <- ncol(X)
          # check if there are enough regressors
          if (m > 1) {
            WX <- matrix(nrow=nrow(X),ncol=(m-(K-1)))
            for (k in K:m) {
              wx <- lag.listw(listw, X[,k], 
                              zero.policy=zero.policy)
              if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
              WX[,(k-(K-1))] <- wx
            }
          }
          if (K == 2) {
            # unnormalized weight matrices
            if (!(listw$style == "W")) {
              intercept <- as.double(rep(1, nrow(X)))
              wx <- lag.listw(listw, intercept, 
                              zero.policy = zero.policy)
              if (m > 1) {
                WX <- cbind(wx, WX)
              } else {
                WX <- matrix(wx, nrow = nrow(X), ncol = 1)
              }
            } 
          }   
          X <- cbind(X, WX)
          #  accommodate aliased coefficients 120314
          if (any(object$aliased))
            X <- X[,-which(object$aliased)]
          trend <- X %*% B
          signal <- rep(0, length(trend)) # vecteur de 0
          res <- trend + signal # prediction
          attr(res, "trend") <- trend
          attr(res, "signal") <- signal
        } else stop("unkown error model etype")
      } else if (object$type == "mixed") { # Wy+WX
        if (is.null(listw) || !inherits(listw, "listw")) 
          stop ("spatial weights list required")
        if (nrow(newdata) != length(listw$neighbours))
          stop("mismatch between newdata and spatial weights")
        B <- object$coefficients
        #			mt <- terms(object$formula, data = newdata)
        frm <- formula(object$call)
        mt <- delete.response(terms(frm, data=newdata))
        #			mf <- lm(object$formula, newdata, method="model.frame")
        mf <- model.frame(mt, newdata)
        X <- model.matrix(mt, mf)
        K <- ifelse(colnames(X)[1] == "(Intercept)", 2, 1)
        m <- ncol(X)
        # check if there are enough regressors
        if (m > 1) {
          WX <- matrix(nrow=nrow(X),ncol=(m-(K-1)))
          for (k in K:m) {
            wx <- lag.listw(listw, X[,k], 
                            zero.policy=zero.policy)
            if (any(is.na(wx))) 
              stop("NAs in lagged independent variable")
            WX[,(k-(K-1))] <- wx
          }
        }
        if (K == 2) {
          # unnormalized weight matrices
          if (!(listw$style == "W")) {
            intercept <- as.double(rep(1, nrow(X)))
            wx <- lag.listw(listw, intercept, 
                            zero.policy = zero.policy)
            if (m > 1) {
              WX <- cbind(wx, WX)
            } else {
              WX <- matrix(wx, nrow = nrow(X), ncol = 1)
            }
          } 
        }   
        X <- cbind(X, WX)
        #  accommodate aliased coefficients 120314
        if (any(object$aliased))
          X <- X[,-which(object$aliased)]
        trend <- X %*% B
        if (power) { # calcul de (I - rho*W)^-1 en élevant à la puissance
          W <- as(listw, "CsparseMatrix")
          res <- c(as(powerWeights(W, rho=object$rho,
                                   X=trend, order=order, tol=tol), "matrix"))
        } else { # calcul de (I - rho*W)^-1 en inversant la matrice
          res <- c(invIrW(listw, object$rho) %*% trend)
        }
        if (legacy) {
          signal <- object$rho * lag.listw(listw, 
                                           res, zero.policy=zero.policy)
          res <- c(trend + signal)
        } else {
          signal <- res - trend
        }
        #                        if (pred.se) {
        #                            samples <- attr(lagImpact, "samples")$samples
        #                            irho <- attr(lagImpact, "samples")$irho
        #                            drop2beta <- attr(lagImpact, "samples")$drop2beta
        #                            nSim <- nrow(samples)
        #                            outmat <- matrix(NA, ncol=nSim, nrow=nrow(X))
        #                            for (i in 1:nSim) {
        #                                B <- samples[i, -drop2beta]
        #                                trend <- X %*% B
        #                                rho <- samples[i, irho]
        #                                if (power) {
        #                                    res <- c(as(powerWeights(W, rho=rho,
        #                                    X=trend, order=order, tol=tol), "matrix"))
        #                                } else {
        #                                    res <- c(invIrW(listw, rho) %*% trend)
        #                                }
        #                                outmat[,i] <- res
        #                            }
        #                            pred.se <- apply(outmat, 1, sd)
        #                            attr(res, "pred.se") <- pred.se
        #                        }
        attr(res, "trend") <- c(trend)
        attr(res, "signal") <- c(signal)
      } else { # Wy (ou Wy+We ??)
        if (is.null(listw) || !inherits(listw, "listw")) 
          stop ("spatial weights list required")
        if (nrow(newdata) != length(listw$neighbours))
          stop("mismatch between newdata and spatial weights")
        B <- object$coefficients
        #			mt <- terms(object$formula, data = newdata)
        frm <- formula(object$call)
        mt <- delete.response(terms(frm, data=newdata))
        #			mt <- delete.response(terms(object$formula))
        #			mf <- lm(object$formula, newdata, method="model.frame")
        # resolved problem of missing response column in newdata reported by
        # Christine N. Meynard, 060201
        mf <- model.frame(mt, newdata)
        if (dim(mf)[1] != length(listw$neighbours))
          stop("missing values in newdata")
        X <- model.matrix(mt, mf)
        #  accommodate aliased coefficients 120314
        if (any(object$aliased))
          X <- X[,-which(object$aliased)]
        trend <- X %*% B
        if (power) {
          W <- as(listw, "CsparseMatrix")
          res <- c(as(powerWeights(W, rho=object$rho,
                                   X=trend, order=order, tol=tol), "matrix"))
        } else {
          res <- c(invIrW(listw, object$rho) %*% trend)
        }
        if (legacy) {
          signal <- object$rho * lag.listw(listw, 
                                           res, zero.policy=zero.policy)
          res <- c(trend + signal)
        } else {
          signal <- res - trend
        }
        #                        if (pred.se) {
        #                            samples <- attr(lagImpact, "samples")$samples
        #                            irho <- attr(lagImpact, "samples")$irho
        #                            drop2beta <- attr(lagImpact, "samples")$drop2beta
        #                            nSim <- nrow(samples)
        #                            outmat <- matrix(NA, ncol=nSim, nrow=nrow(X))
        #                            for (i in 1:nSim) {
        #                                B <- samples[i, -drop2beta]
        #                                trend <- X %*% B
        #                                rho <- samples[i, irho]
        #                                if (power) {
        #                                    res <- c(as(powerWeights(W, rho=rho,
        #                                    X=trend, order=order, tol=tol), "matrix"))
        #                                } else {
        #                                    res <- c(invIrW(listw, rho) %*% trend)
        #                                }
        #                                outmat[,i] <- res
        #                            }
        #                            pred.se <- apply(outmat, 1, sd)
        #                            attr(res, "pred.se") <- pred.se
        #                        }
        attr(res, "trend") <- c(trend)
        attr(res, "signal") <- c(signal)
      }
    } else { 
      if (type %in% c("TS1", "KP2")) { # TODO: add "X" type for out-of-sample. Becareful with mixed models: Wss %*% X is include in object$X
        if (nrow(newdata) != 1) # only in leave-one-out
          stop("predictor type only for leave-one-out (newdata should have one row)") # TODO: should we iterate them?
        is.newdata <- attr(listw$neighbours, "region.id") == rownames(newdata)
        WoS <- as(listw, "CsparseMatrix")[ is.newdata, ! is.newdata ] # TODO: with listw$neighbours and listw$weights
        TS1 <- newdata %*% B + object$rho + WoS %*% y
        if (type == "TS1") {
          res <- TS1
          #TODO: trend and signal
        } else { # KP2
          #TODO
          #KP2 <- TS1 + 
        }
      } else { # not leave-one-out
        listw.d = .listw.decompose(listw, region.id.data, region.id.newdata, type = c("Wss", "Wos", "Wso", "Woo"))
        
      }
    }
  }
  class(res) <- "sarlm.pred"
  res
}

# decompose a listw object into Wss Wso Wos and Woo sparse matrices
.listw.decompose <- function(listw, region.id.data, region.id.newdata, type = c("Wss", "Wos", "Wso", "Woo")) { # TODO: hidden? in this file? zero.policy?
  if (is.null(listw) || !inherits(listw, "listw")) 
    stop ("spatial weights list required")
  region.id <- attr(listw, "region.id")
  if (!all(region.id.data %in% region.id))
    stop("at least one region.id in data is not in listw object")
  if (!all(region.id.newdata %in% region.id))
    stop("at least one region.id in newdata is not in listw object")
  if (!all(type %in% c("Wss", "Wos", "Wso", "Woo"))) #TODO: better way to do?
    stop("type is incorrect")
  W <- as(listw, "CsparseMatrix")
  s <- list(Wss = NULL, Wos = NULL, Wso = NULL, Woo = NULL)
  if ("Wss" %in% type)
    s$Wss <- W[region.id.data, region.id.data]
  if ("Wos" %in% type)
    s$Wos <- W[region.id.newdata, region.id.data]
  if ("Wso" %in% type)
    s$Wso <- W[region.id.data, region.id.newdata]
  if ("Woo" %in% type)
    s$Woo <- W[region.id.newdata, region.id.newdata]
  return(s)
  # TODO: be careful when region.id names cannot be the names of a sparse matrix
  # becareful to the order of rows and columns in the sparse matrix when we used it
}
print.sarlm.pred <- function(x, ...) {
	res <- as.data.frame(x)
	print(res, ...)
	invisible(res)
}


as.data.frame.sarlm.pred <- function(x, ...) {
    res <- data.frame(fit=as.vector(x), trend=attr(x, "trend"), 
        signal=attr(x, "signal"))
    res
}

