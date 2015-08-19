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


# KP2 / KP3
if (object$type %in% c("sac", "sacmixed")) {
  if (power){
    GL <- powerWeights(Wi, rho= object$lambda, order= order, tol= tol, X= Diagonal(length(ys)+1))
    GR <- powerWeights(Wi, rho= object$rho, order= order, tol= tol, X= Diagonal(length(ys)+1))
  } else {
    GL <- invIrW(Wi, object$lambda)
    GR <- invIrW(Wi, object$rho)
  }
  sum.u <- GL %*% t(GL)
  sum.y <- GR %*% sum.u %*% t(GR)
} else if (object$type %in% c("lag", "mixed")) {
  if (power){
    GR <- powerWeights(Wi, rho= object$rho, order= order, tol= tol, X= Diagonal(length(ys)+1))
  } else {
    GR <- invIrW(Wi, object$rho)
  }
  sum.u <- Diagonal(length(ys)+1)
  sum.y <- GR %*% t(GR)
} else if (object$type == "error") {
  if (power){
    GL <- powerWeights(Wi, rho= object$lambda, order= order, tol= tol, X= Diagonal(length(ys)+1))
  } else {
    GL <- invIrW(Wi, object$lambda)
  }
  sum.u <- GL %*% t(GL)
  sum.y <- sum.u
} else stop("unknown model type")


if (type == "KP2") {
  region.id.data <- attr(ys, "names")
  region.id.newdata <- rownames(newdata)
  W <- as(listw, "CsparseMatrix")
  KP2 <- rep(NA, nrow(newdata))
  for (i in 1:nrow(newdata)) {
    region.id.temp <- c(region.id.data, region.id.newdata[i])
    Wi <- W[region.id.temp, region.id.temp]
    Xi <- rbind(Xs, Xo[i,])
    wi <- Wi[length(region.id.temp), ]
    yi <- c(ys, 0)
    covar <- sum.u[length(region.id.temp), ] %*% t(GR) %*% t(wi) / (wi %*% sum.y %*% t(wi))
    Ewiy <- wi %*% GR %*% Xi %*% B
    KP2[i] <- TS1[i] + covar %*% (wi %*% yi - Ewiy)
  }
}

if (type == "KP3") {
  region.id.data <- attr(ys, "names")
  region.id.newdata <- rownames(newdata)
  W <- as(listw, "CsparseMatrix")
  KP2 <- rep(NA, nrow(newdata))
  for (i in 1:nrow(newdata)) {
    region.id.temp <- c(region.id.data, region.id.newdata[i])
    Wi <- W[region.id.temp, region.id.temp]
    Xi <- rbind(Xs, Xo[i,])
    wi <- Wi[length(region.id.temp), ]
    yi <- c(ys, 0)
    cov <- sum.u[length(region.id.temp), ] %*% t(GR)[, -length(region.id.temp)]
    KP3[i] <- TS1[i] + cov %*% solve(sum.y[-length(region.id.temp), -length(region.id.temp)]) %*% (ys - GR[-length(region.id.temp),] %*% Xi %*% B)
  }
}

# KP5

