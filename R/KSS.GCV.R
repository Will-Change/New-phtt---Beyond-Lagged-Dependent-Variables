FUN.iterate.GCV <- function(TR.Y.mat, TR.X.mat, N, T, P){
	y     <- matrix(TR.Y.mat,  nrow = N*T, ncol = 1)                  # (TN x 1)
	x     <- matrix(TR.X.mat,  nrow = N*T, ncol = P)                  # (TN x P)
	t.seq <- seq(0, 1, length.out = T)

	xx         <- crossprod(x)                               # (p x p)
	L <- max(eigen(xx)$values) / N * T                        # Lipschitz constant
	inv.xx     <- solve(xx)                                  # (p x p)
	inv.xx.x   <- inv.xx %*% t(x)                            # (p x TN)

	# Define the functions as before
	FUN.ols.beta <- function(updated.y, x, inv.xx.x){
		beta <- tcrossprod(inv.xx.x, t(updated.y))
	}

	FUN.lasso.beta <- function(updated.y, x) {
		# Fit Lasso model; note glmnet takes y as a vector and x as a matrix
		lasso.model <- glmnet(x, updated.y, alpha = 1)

		# Use cross-validation to find optimal lambda
		cv.lasso <- cv.glmnet(x, updated.y, alpha = 1)
		best.lambda <- cv.lasso$lambda.min

		# Extract coefficients at best lambda
		beta <- as.matrix(coef(lasso.model, s = best.lambda)[-1, , drop = FALSE]) # Dropping intercept
		return(beta)
	}

	# Starting values for FISTA
	ymats <- matrix(y, T, N)
	xmats <- matrix(x, T, (N*P))
	zmats <- cbind(ymats, xmats)
	trzma <- zmats - svd.pca(zmats, given.d = round(sqrt(min(N, T))))$Q.fit
	tryxm <- matrix(trzma, (T*N), (P+1))
	try   <- tryxm[, 1, drop = FALSE]
	trx   <- tryxm[, -1, drop = FALSE]
	beta.0 <- coef(lm(try ~ -1 + trx))                      # (p x p)

	t_k <- 1
	y_k <- beta.0  # Initialize y_k with beta_0

	## FISTA Iteration
	inner.iteration <- function(y, x, inv.xx.x = inv.xx.x, beta.0 = beta.0, t.seq, t_k, y_k, i = 1){
		nr <- length(t.seq)
		nc <- nrow(y) / nr

		# Update using y_k
		w.0 <- y - tcrossprod(x, t(beta.0))
		W.0 <- matrix(w.0, nr, nc)
		PCA.0 <- smooth.Pspline(x = t.seq, y = W.0, method = 3)
		y.fitted.0 <- PCA.0$ysmth
		y.updated.1 <- y - c(y.fitted.0)
		beta.1 <- FUN.lasso.beta(y.updated.1, x)

		# FISTA updates
		t_k_plus_1 <- (1 + sqrt(1 + 4 * t_k^2)) / 2
		y_k_plus_1 <- beta.1 + ((t_k - 1) / t_k_plus_1) * (beta.1 - beta.0)

		# Convergence condition
		if (all(abs(beta.0 - beta.1) < 1e-3) | i == 100) {
			Result <- list(PCA = PCA.0, beta = beta.1, Nbr.Iterations = i)
			return(Result)
		} else {
			return(inner.iteration(y = y, x = x, inv.xx.x = inv.xx.x, beta.0 = beta.1, t.seq = t.seq, t_k = t_k_plus_1, y_k = y_k_plus_1, i = i + 1))
		}
	}

	Result <- inner.iteration(y = y, x = x, t.seq = t.seq, t_k = t_k, y_k = y_k, beta.0 = beta.0, i = 1)
	return(Result)
}
