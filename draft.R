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
TS <- trend + object$rho * W %*% y


# out-of-sample


