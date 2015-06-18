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
predict.sarlm <- function(object, newdata=NULL, listw=NULL, type=NULL, all.data=FALSE,
                          zero.policy=NULL, legacy=TRUE, power=NULL, order=250, tol=.Machine$double.eps^(3/5), #pred.se=FALSE, lagImpact=NULL, 
                          ...) {
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (object$type == "sac") stop("no predict method for sac") # TODO: need to change # Wy et We
  if (is.null(power)) power <- object$method != "eigen"
  stopifnot(is.logical(all.data))
  stopifnot(is.logical(legacy))
  stopifnot(is.logical(power))
  #        if (pred.se && object$type == "error") {
  #            pred.se <- FALSE
  #            warning("standard error estimates not available for error models")
  #        }
  #        if (pred.se && is.null(lagImpact))
  #            stop("lagImpact object from impact method required for standard error estimate")
  Xs <- object$X
  B <- object$coefficients
  ys <- object$y
  tarXs <- object$tarX
  tarys <- object$tary
  trends <- Xs %*% B
  
  if (is.null(newdata)) { # in-sample pred
    if (is.null(type) || (type == "TS")) { # defaut predictor
      res <- fitted.values(object)
      if (object$type == "error") { # We ou WX+We
        attr(res, "trend") <- as.vector(trends)
        attr(res, "signal") <- as.vector( -1 * (tarys - ys) -	-1 * (tarXs - Xs) %*% B)
      } else { # lag model or other
        attr(res, "trend") <- as.vector(trends)
        attr(res, "signal") <- as.vector( -1 * (tarys - ys))
      }
    } else { # new predictors
      if (! type %in% c("X", "TC", "BP")) stop("no such predictor type")
      #TODO: lag model only? + adapt to SDM model
      #TODO: match newdata columns order with formula
      if (is.null(listw) || !inherits(listw, "listw"))
        stop ("spatial weights list required")
      if (nrow(Xs) != length(listw$neighbours))
        stop("mismatch between data and spatial weights")
      if (type == "X") { # warning: equals to WX predictor for mixed model
        res <- as.vector(trends)
      } else if (type %in% c("TC", "BP")) { # need to compute TC
        if (power){
          W <- as(listw, "CsparseMatrix")
          TC <- powerWeights(W, rho = object$rho, X = Xs, order = order, tol = tol) %*% B # TODO: not good for SDM model!
        } else {
          TC <- invIrW(listw, object$rho) %*% trends # TODO: not good for SDM model!
        }
      }
      if(type == "TC") res <- as.vector(TC)
      if(type == "BP") {
        W <- as(listw, "CsparseMatrix")
        Qss <- 1/object$s2 * (Diagonal(dim(W)[1]) - (object$rho * t(W))) %*% (Diagonal(dim(W)[1]) - (object$rho * W)) # precision matrix for LAG model
        DiagQss <- Diagonal(x = diag(Qss))
        BP <- TC - solve(DiagQss) %*% (Qss - DiagQss) %*% (ys - TC)
        # Can BP also be applied to the SEM model? Cf LeSage and Pace (2004). Note: \hat{\mu_i} need to be adapted
        res <- as.vector(BP)
      }
      attr(res, "trend") <- as.vector(trends) # TODO: to be removed?
      attr(res, "signal") <- NULL
    }
  }
  else { # out-of-sample
    #CHECK
    if (!((is.null(type) || (type == "TS")) && object$type == "error" && object$etype == "error")) { # need of listw (ie. not in the case of defaut predictor and SEM model)
      if (is.null(listw) || !inherits(listw, "listw"))
        stop ("spatial weights list required")
      if (any(! rownames(newdata) %in% attr(listw, "region.id")))
        stop("mismatch between newdata and spatial weights")
      if ((is.null(type) || (type == "TS"))) { # only need Woo
        if (any(! attr(listw, "region.id") %in% rownames(newdata))) { # for consistency, allow the use of a larger listw with in-sample data
          listw <- subset.listw(listw, (attr(listw, "region.id") %in% rownames(newdata)), zero.policy = zero.policy) # TODO: need more test
          # TODO: reorder listw if newdata is not in the order of listw
        }
      }
      #...
    }
    
    # DATA
    frm <- formula(object$call)
    mt <- delete.response(terms(frm, data=newdata)) # returns a terms object for the same model but with no response variable
    mf <- model.frame(mt, newdata)
    # resolved problem of missing response column in newdata reported by
    # Christine N. Meynard, 060201
    if (dim(mf)[1] != length(listw$neighbours))
      stop("missing values in newdata")
    Xo <- model.matrix(mt, mf)
    
    if (object$type == "mixed" || object$etype == "mixed" ) { # mixed model: compute WXo
      K <- ifelse(colnames(Xo)[1] == "(Intercept)", 2, 1)
      m <- ncol(Xo)
      # check if there are enough regressors
      if (m > 1) {
        WXo <- matrix(nrow=nrow(Xo),ncol=(m-(K-1)))
        for (k in K:m) {
          wx <- lag.listw(listw, Xo[,k], 
                          zero.policy=zero.policy)
          if (any(is.na(wx)))
            stop("NAs in lagged independent variable")
          WXo[,(k-(K-1))] <- wx
        }
      }
      if (K == 2) {
        # unnormalized weight matrices
        if (!(listw$style == "W")) {
          intercept <- as.double(rep(1, nrow(Xo)))
          wx <- lag.listw(listw, intercept, 
                          zero.policy = zero.policy)
          if (m > 1) {
            WXo <- cbind(wx, WXo)
          } else {
            WXo <- matrix(wx, nrow = nrow(Xo), ncol = 1)
          }
        } 
      }   
      Xo <- cbind(Xo, WXo)
      #  accommodate aliased coefficients 120314
      if (any(object$aliased))
        Xo <- Xo[,-which(object$aliased)]
      trendo <- Xo %*% B
    } else {
      #  accommodate aliased coefficients 120314
      if (any(object$aliased))
        Xo <- Xo[,-which(object$aliased)]
      trendo <- Xo %*% B
    }
    
    if (is.null(type) || (type == "TS")) { # defaut predictor
      if (object$type == "error") {
        if (object$etype == "error") { # We
          signal <- rep(0, length(trendo))
          res <- trendo + signal
          attr(res, "trend") <- trendo
          attr(res, "signal") <- signal
        } else if (object$etype == "emixed") { # WX + We
          signal <- rep(0, length(trendo))
          res <- trend + signal
          attr(res, "trend") <- trendo
          attr(res, "signal") <- signal
        } else stop("unkown error model etype")
      } else if (object$type == "mixed") { # Wy+WX
        if (power) {
          W <- as(listw, "CsparseMatrix")
          res <- c(as(powerWeights(W, rho=object$rho,
                                   X=trendo, order=order, tol=tol), "matrix"))
        } else { # calcul de (I - rho*W)^-1 en inversant la matrice
          res <- c(invIrW(listw, object$rho) %*% trendo)
        }
        if (legacy) {
          signal <- object$rho * lag.listw(listw, 
                                           res, zero.policy=zero.policy)
          res <- c(trendo + signal)
        } else {
          signal <- res - trendo
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
        attr(res, "trend") <- c(trendo)
        attr(res, "signal") <- c(signal)
      } else { # Wy (TODO: when we will allow predict for SAC, make sure we have an error here)
        if (power) {
          W <- as(listw, "CsparseMatrix")
          res <- c(as(powerWeights(W, rho=object$rho,
                                   X=trendo, order=order, tol=tol), "matrix"))
        } else {
          res <- c(invIrW(listw, object$rho) %*% trendo)
        }
        if (legacy) {
          signal <- object$rho * lag.listw(listw, 
                                           res, zero.policy=zero.policy)
          res <- c(trendo + signal)
        } else {
          signal <- res - trendo
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
        attr(res, "trend") <- c(trendo)
        attr(res, "signal") <- c(signal)
      }
    } else {
      if (type %in% c("TS1", "KP2")) {
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        Wos <- .listw.decompose(listw, region.id.data = attr(y, "names"), region.id.newdata = rownames(newdata), type = "Wos")$Wos
        TS1 <- Xo %*% B + object$rho * Wos %*% ys # TODO: what to todo when a same spatial unit is on newdata and on data
        if (type == "TS1") {
          res <- TS1
          #TODO: trend and signal
        } else { # KP2
          #TODO
          #KP2 <- TS1 + 
        }
      } else { # not leave-one-out
        if (type == "TC") {
          #notations of C.Thomas and al (2015)
          if (all.data) { # TCo = TC for out-of-sample spatial units
          listw.d = .listw.decompose(listw, region.id.data = NULL, region.id.newdata = rownames(newdata), type = c("Wss", "Wos", "Wso", "Woo"))
          Wss <- listw.d$Wss
          Wso <- listw.d$Wso
          Wos <- listw.d$Wos
          Wss <- listw.d$Wss
          mB <- - object$rho * Wso
          mC <- - object$rho * Wos
          mD <- Diagonal(dim(Woo)[1]) - object$rho * Woo
          if (power){
            mAInvXsB <- powerWeights(Wss, rho = object$rho, X = Xs, order = order, tol = tol) %*% B
            mAInvmB <- powerWeights(Wss, rho = object$rho, X = mB, order = order, tol = tol)
            E <- solve(mD - mC %*% mAInvmB)
            TCo <- - E %*% mC %*% mAInvXsB + E %*% Xo %*% B
          } else {
            # manually without using invIrW, because it use method="solve" by default, but accept only a listw object. 
            mA <- Diagonal(dim(Wss)[1]) - object$rho * Wss
            mAInv <- solve(mA)
            E <- solve(mD - mC %*% mAInv %*% mB)
            TCo <- - E %*% mC %*% mAInv %*% Xs %*% B + E %*% Xo %*% B
          }
          res <- TCo
          } else {
            # compute s and o units together
            X <- rbind(Xs, Xo)
            trend <- X %*% B
            if (power){
              W <- as(listw, "CsparseMatrix")
              TC <- powerWeights(W, rho = object$rho, X = X, order = order, tol = tol) %*% B
            } else {
              TC <- invIrW(listw, object$rho) %*% trend
            }
            res <- TC
            # TODO: performances test to know if computing TCo and TCs is quicker than computing TC
          }
          
        } else {
          stop("unknow predictor type")
        }
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
  # TODO: what to do if region.id.newdata is included in region.id.data. Include it, or remove it from data?
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

