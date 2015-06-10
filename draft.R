# implementation of 'About predictions in spatial autoregressive models: optimal and almost optimal strategies'

# in-sample
if (power){
   W <- as(listw, "CsparseMatrix")
   TC <- powerWeights(W, rho = rho, X = X, order = order, tol = tol) %*% B
} else {
   TC <- invIrW(listw, rho) %*% X %*% B
}

W <- as(listw, "CsparseMatrix")
Qs <- 1/object$s2 * (as(diag(dim(W)[1]), "CsparseMatrix") - (rho * t(W))) %*% (as(diag(dim(W)[1]), "CsparseMatrix") - (rho * W)) # precision matrix for LAG model # TODO : change to more appropriate type of Matrix: diagonalMatrix?
BP <- TC - solve(as(diag(diag(Qs)), "CsparseMatrix")) %*% (Qs - diag(Qs)) %*% (y - TC) # TODO : change to more appropriate type of Matrix: diagonalMatrix?

W <- as(listw, "CsparseMatrix")
TS <- X %*% B + rho * W %*% y

# out-of-sample


