# no need:
#X or WX -> trend

#needed:

#KP1 -> (I-rho*W)^-1 %*% Xb -> TC1
#TODO: use the definition of TCo of C.Thomas applied in the leave-one-out case: quicker?
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
    n <- dim(Wi)[1]
    mat <- diag(n) - object$rho * Wi
    trendi <- c(trends, trendo[i])
    res[i] <- (solve(mat) %*% trendi)[length(region.id.temp)]
  }
}

#XWy
#TODO: Where it comes from???
#TODO: new name?
if (object$type == "error") {
  prdWXy <- trends + object$lambda* lag.listw(listw, ys - trends)
} else if (object$type %in% c("lag", "mixed")) {
  prdWXy <- trends + object$rho* lag.listw(listw, ys)
} else if (object$type %in% c("sac", "sacmixed")) {
  prdWXy <- trends + object$rho* lag.listw(listw, ys) + object$lambda* lag.listw(listw, ys - trends)
} else {
  stop("unkown model type")
}


prdWXy <- trends + object$rho* lag.listw(listw, ys) + object$lambda* lag.listw(listw, ys - trends)

# LSP


# KPG


# KP2


# KP3


