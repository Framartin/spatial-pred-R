# implementation of 'About predictions in spatial autoregressive models: optimal and almost optimal strategies'

# TODO; optimize using tcrossprod() when possible

# in-sample
if (power){
   W <- as(listw, "CsparseMatrix")
   TC <- powerWeights(W, rho = object$rho, X = X, order = order, tol = tol) %*% B
} else {
   TC <- invIrW(listw, object$rho) %*% trend #X %*% B
}

W <- as(listw, "CsparseMatrix")
Qss <- 1/object$s2 * (as(diag(dim(W)[1]), "CsparseMatrix") - (object$rho * t(W))) %*% (as(diag(dim(W)[1]), "CsparseMatrix") - (object$rho * W)) # precision matrix for LAG model # TODO : change to more appropriate type of Matrix: diagonalMatrix?
DiagQss = as(diag(diag(Qss)), "CsparseMatrix") # TODO : change to more appropriate type of Matrix: diagonalMatrix?
BP <- TC - solve(DiagQss) %*% (Qss - DiagQss) %*% (y - TC)
# Can BP also be applied to the SEM model? Cf LeSage and Pace (2004). Note: \hat{\mu_i} need to be adapted

W <- as(listw, "CsparseMatrix")
TS <- trend + object$rho * W %*% y # correspond to the object$fitted.values for a lag model


# out-of-sample

## Leave-one-out
TS1 <- 

## other cases

# TC

#mA <- Diagonal(dim(Wss)[1]) - object$rho * Wss
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
# TODO: TCs

# compute s and o units together
X <- rbind(Xs, Xo)
trend <- X %*% B
if (power){
   W <- as(listw, "CsparseMatrix")
   TC <- powerWeights(W, rho = object$rho, X = X, order = order, tol = tol) %*% B
} else {
   TC <- invIrW(listw, object$rho) %*% trend
}
# TODO: performances test


# BP
TCo <- TC[is.newdata]
TCs <- TC[is.data]
Q <- 1/object$s2 * ( Diagonal(dim(W)[1]) - object$rho * (t(W) + W) + object$rho^2 * (t(W) %*% W) )
Qoo <- Q[is.newdata, is.newdata]
Qos <- Q[is.newdata, is.data]
BPo <- TCo - solve(Qoo) %*% Qos %*% (ys - TCs)

#BPw

if (power){
   W <- as(listw, "CsparseMatrix")
   invW <- powerWeights(W, rho = object$rho, X = Diagonal(dim(W)[1]), order = order, tol = tol)
} else {
   invW <- invIrW(listw, object$rho)
}
X <- rbind(Xs, Xo)
TC <- invW %*% X %*% B
is.data <- 1:length(ys)
is.newdata <- (length(ys)+1):length(TC)
TCo <- TC[is.newdata]
TCs <- TC[is.data]
#Sigma <- object$s2 * solve((Diagonal(dim(W)[1]) - object$rho * t(W)) %*% (Diagonal(dim(W)[1]) - object$rho * W))
Sigma <- object$s2 * invW %*% t(invW)
Sos <- Sigma[is.newdata, is.data]
Sss <- Sigma[is.data, is.data]
Wos <- .listw.decompose(listw, region.id.data = attr(ys, "names"), region.id.newdata = rownames(newdata), type = "Wos")$Wos
BPW <- TCo + Sos %*% t(Wos) %*% solve(Wos %*% Sss %*% t(Wos)) %*% (Wos %*% ys - Wos %*% TCs)


# BPn
#TODO: approx of BP quicker to compute: replace S by J, the spatial units of S which are neighbourgs of O
is.data <- 1:length(ys)
is.newdata <- (length(ys)+1):length(TC)
TCs <- TC[is.data]
TCo <- TC[is.newdata]
# compute J = set of all sites in S which are neighbors of at least one site in O
O <- which(attr(listw,"region.id") %in% rownames(newdata))
S <- which(attr(listw,"region.id") %in% attr(ys, "names"))
J.logical <- rep(FALSE, length(listw$neighbours))
for (i in S) {
  J.logical[i] <- any(O %in% listw$neighbours[[i]])
}
J <- attr(listw,"region.id")[J.logical]

if (length(J)<1) {
  warning("out-of-sample units have no neighbours")
  BPN <- TCo
} else {
  W <- as(listw, "CsparseMatrix")
  region.id <- c(J, rownames(newdata))
  W_jo <- W[region.id, region.id]
  rm(W)
  Q_jo <- 1/objects$s2 * (Diagonal(length(region.id)) - object$rho * (W_jo + t(W_jo)) + object$rho^2 * (t(W_jo) %*% W_jo))
  is.j <- 1:length(J)
  is.o <- (length(J)+1):length(region.id)
  Qoo <- Q_jo[is.o, is.o]
  Qoj <- Q_jo[is.o, is.j]
  rm(Q_jo)
  yj <- ys[J]
  TCj <- TCs[attr(ys, "names") %in% J]
  BPN <- as.vector(TCo - solve(Qoo) %*% Qoj %*% (yj - TCj))
}

