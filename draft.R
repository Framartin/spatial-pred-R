# implementation of 'About predictions in spatial autoregressive models: optimal and almost optimal strategies'

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

