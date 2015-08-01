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
                          spChk=NULL, ...) {
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (is.null(type)) type <- "default"
  # check type with model
  if (type %in% c("default", "TS") & object$type %in% c("sac", "sacmixed")) stop("no such predict method for sac model")
  if (type %in% c("TC", "TS", "BP", "TS1") & object$type == "error") stop("no such predict method for error model")
  if (type %in% c("TC", "TS", "BP", "TS1") & object$type %in% c("sac", "sacmixed")) warning("predict method developed for lag model, use carefully")
  
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
    if (type == "default" || type == "TS") { # defaut predictor
      res <- fitted.values(object)
      if (object$type == "error") { # We ou WX+We
        attr(res, "trend") <- as.vector(trends)
        attr(res, "signal") <- as.vector( -1 * (tarys - ys) -	-1 * (tarXs - Xs) %*% B)
      } else { # lag model or other
        attr(res, "trend") <- as.vector(trends)
        attr(res, "signal") <- as.vector( -1 * (tarys - ys))
      }
    } else { # new predictors
      if (type != "trend") { # need listw
        if (is.null(listw) || !inherits(listw, "listw"))
          stop ("spatial weights list required")
        if (nrow(Xs) != length(listw$neighbours))
          stop("mismatch between data and spatial weights")
        if (is.null(spChk)) spChk <- get.spChkOption()
        if (spChk && !chkIDs(Xs, listw))
          stop("Check of data and weights ID integrity failed")
      }
      
      if (type %in% c("TC", "BP")) { # need to compute TC
        if (power){
          W <- as(listw, "CsparseMatrix")
          TC <- powerWeights(W, rho = object$rho, X = Xs, order = order, tol = tol) %*% B
        } else {
          TC <- invIrW(listw, object$rho) %*% trends
        }
      }
      
      if (type == "trend") {
        res <- as.vector(trends)
      } else if(type == "TC") {
        res <- as.vector(TC)
      } else if(type == "BP") {
        W <- as(listw, "CsparseMatrix")
        Qss <- 1/object$s2 * (Diagonal(dim(W)[1]) - (object$rho * t(W))) %*% (Diagonal(dim(W)[1]) - (object$rho * W)) # precision matrix for LAG model
        DiagQss <- Diagonal(x = diag(Qss))
        BP <- TC - solve(DiagQss) %*% (Qss - DiagQss) %*% (ys - TC)
        # TODO: Can BP also be applied to the SEM model? Cf LeSage and Pace (2004). Note: \hat{\mu_i} need to be adapted
        res <- as.vector(BP)
      } else {
        stop("no such in-sample predictor type")
      }
      attr(res, "trend") <- as.vector(trends) # TODO: to be removed?
      #attr(res, "signal") <- NULL
    }
  } else { # out-of-sample
    #CHECK
    if (is.null(rownames(newdata))) stop("newdata should have region.id as rownames")
    if (!(type == "default" && object$type == "error" && object$etype == "error") && type != "trend") { # need of listw (ie. neither in the case of defaut predictor and SEM model, nor trend type)
      if (is.null(listw) || !inherits(listw, "listw"))
        stop ("spatial weights list required")
      if (any(! rownames(newdata) %in% attr(listw, "region.id")))
        stop("mismatch between newdata and spatial weights. newdata should have region.id as rownames")
      
      #TODO: add the case when we decompose W in Woo, Wso, Wos and Wss
      
      # wanted order of the listw
      if (type == "default") { # only need Woo
        region.id <- rownames(newdata)
      } else {
        region.id <- c(attr(ys, "names"), rownames(newdata))
      }
      
      if (length(region.id) != length(attr(listw, "region.id")) || !all(region.id == attr(listw, "region.id"))) { # if listw is not directly ok
        if (all(subset(attr(listw, "region.id"), attr(listw, "region.id") %in% region.id) == region.id)) { # only need a subset.listw, ie. spatial units are in the right order
          listw <- subset.listw(listw, (attr(listw, "region.id") %in% rownames(newdata)), zero.policy = zero.policy)
          W <- as(listw, "CsparseMatrix")
        } else { # we use a sparse matrix transformation to reorder a listw
          W <- as(listw, "CsparseMatrix")
          W <- W[region.id, region.id]
          style <- listw$style
          listw <- mat2listw(W, row.names = region.id, style = style) # re-normalize to keep the style
          W <- as(listw, "CsparseMatrix")
        }
      }
      #optional check
      if (is.null(spChk)) spChk <- get.spChkOption()
      names(region.id) <- region.id
      if (spChk && !chkIDs(temp, listw))
        stop("Check of data and weights ID integrity failed")
    }
    

    
    # DATA
    frm <- formula(object$call)
    mt <- delete.response(terms(frm, data=newdata)) # returns a terms object for the same model but with no response variable
    mf <- model.frame(mt, newdata)
    # resolved problem of missing response column in newdata reported by
    # Christine N. Meynard, 060201
    if (dim(mf)[1] != nrow(newdata))
      stop("missing values in newdata")
    Xo <- model.matrix(mt, mf)
    
    if (sum(object$type == "mixed", object$etype == "mixed" )>0) { # mixed model: compute WXo # Fix bug if object$etype doesn't exist
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
    
    if (type == "default") { # defaut predictor # TODO: WARNING this code used Woo != C.Thomas TS1 predictor !
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
      } else { # Wy
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
    } else { # new predictors
      if (type %in% c("TS1", "KP4", "KP2")) { # need to compute TS1/KP4
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        Wos <- .listw.decompose(listw, region.id.data = attr(ys, "names"), region.id.newdata = rownames(newdata), type = "Wos")$Wos
        TS1 <- Xo %*% B + object$rho * Wos %*% ys # TODO: what to todo when a same spatial unit is on newdata and on data
      }
      if (type %in% c("TC", "BP")) { # need to compute TC
        #notations of C.Thomas and al (2015)
        if (all.data | type == "BP") {
          # compute s and o units together
          X <- rbind(Xs, Xo)
          trend <- X %*% B
          if (power){
            W <- as(listw, "CsparseMatrix")
            TC <- powerWeights(W, rho = object$rho, X = X, order = order, tol = tol) %*% B
          } else {
            TC <- invIrW(listw, object$rho) %*% trend
          }
          # TODO: performances test to know if computing TCo and TCs is quicker than computing TC
        } else { # TCo = TC for out-of-sample spatial units
          listw.d = .listw.decompose(listw, region.id.data = attr(ys, "names"), region.id.newdata = rownames(newdata), type = c("Wss", "Wos", "Wso", "Woo"))
          Wss <- listw.d$Wss
          Wso <- listw.d$Wso
          Wos <- listw.d$Wos
          Woo <- listw.d$Woo
          mB <- - object$rho * Wso
          mC <- - object$rho * Wos
          mD <- Diagonal(dim(Woo)[1]) - object$rho * Woo
          if (power){
            mAInvXsB <- powerWeights(Wss, rho = object$rho, X = Xs, order = order, tol = tol) %*% B
            mAInvmB <- powerWeights(Wss, rho = object$rho, X = mB, order = order, tol = tol)
            E <- solve(mD - mC %*% mAInvmB)
            TCo <- - E %*% mC %*% mAInvXsB + E %*% Xo %*% B
          } else {
            #mA <- Diagonal(dim(Wss)[1]) - object$rho * Wss
            #mAInv <- solve(mA)
            mAInv <- invIrW(Wss, object$rho)
            E <- solve(mD - mC %*% mAInv %*% mB)
            TCo <- - E %*% mC %*% mAInv %*% Xs %*% B + E %*% Xo %*% B
          }
        }
      }
      
      if (type == "trend") {
        res <- as.vector(trendo)
      } else if (type %in% c("TS1", "KP4")) {
        res <- as.vector(TS1)
      } else if (type == "TC") {
        if(all.data) {
          res <- as.vector(TC)
        } else {
          res <- as.vector(TCo)
        }
      } else if (type == "BP") {
        is.data <- 1:length(ys)
        is.newdata <- (length(ys)+1):length(TC)
        TCo <- TC[is.newdata] # TODO: need to check after modifs on TC: in-sample should be first and then out-of-sample 
        TCs <- TC[is.data] # TODO: idem
        W <- as(listw, "CsparseMatrix")
        Q <- 1/object$s2 * ( Diagonal(dim(W)[1]) - object$rho * (t(W) + W) + object$rho^2 * (t(W) %*% W) )
        Qoo <- Q[is.newdata, is.newdata]
        Qos <- Q[is.newdata, is.data]
        BPo <- TCo - solve(Qoo) %*% Qos %*% (ys - TCs)
        res <- as.vector(BPo)
      } else if (type %in% c("TC1", "KP1")) {
        #TODO: use the definition of TCo of C.Thomas applied in the leave-one-out case: quicker?
        if (nrow(newdata) > 1)
          warning("newdata have more than 1 row and the predictor type is leave-one-out")
        region.id.data <- attr(ys, "names")
        region.id.newdata <- rownames(newdata)
        res <- rep(NA, nrow(newdata))
        W <- as(listw, "CsparseMatrix")
        for (i in 1:nrow(newdata)) {
          region.id.temp <- c(region.id.data, region.id.newdata[i])
          Wi <- W[region.id.temp, region.id.temp]
          Xi <- rbind(Xs, Xo[i,])
          if (power) {
            res[i] <- c(as(powerWeights(Wi, rho=object$rho, X=Xi, order=order, tol=tol), "matrix") %*% B)[length(region.id.temp)]
          } else {
            trendi <- c(trends, trendo[i])
            res[i] <- (invIrW(Wi, object$rho) %*% trendi)[length(region.id.temp)]
          }
        }
      } else {
        stop("unknow predictor type")
      }
    }
  }
  attr(res, "type") <- type
  attr(res, "call") <- match.call()
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
  if (!all(type %in% c("Wss", "Wos", "Wso", "Woo")))
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
  #    res <- data.frame(fit=as.vector(x), trend=attr(x, "trend"), 
  #        signal=attr(x, "signal"))
  #fix bug when no signal or trend attributes
  res <- data.frame(fit=as.vector(x))
  if(!is.null(attr(x, "trend"))) res$trend <- attr(x, "trend")
  if(!is.null(attr(x, "signal"))) res$signal <- attr(x, "signal")
  res
}

