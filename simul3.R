

# Set parameters
n <- 300  # number of individuals
T <- 30  # number of time periods
p <- 2   # number of regressors

# Initialize matrices
X <- array(0, dim = c(T, n, p))
V <- matrix(0, nrow = n, ncol = T)

# Parameter values
R <- matrix(c(0.4, 0.05, 0.05, 0.4), nrow = 2)
mu <- list(c(5, 5), c(7.5, 7.5), c(10, 10))

# Generate Xit values using VAR
for (i in 1:n) {
	  Xi <- as.numeric(solve(diag(2) - R %*% t(R)) %*% matrix(rnorm(2), nrow = 2))
  X[1, i, ] <- Xi + mu[[i %% 3 + 1]]
    for (t in 2:T) {
	        eta <- matrix(rnorm(2), nrow = 2)
      X[t, i, ] <- R %*% matrix(X[t - 1, i, ], nrow = 2) + eta
        }
}

# Generate time-varying individual effects vi(t) for DGP1
for (i in 1:n) {
	  theta <- rnorm(3)
  for (t in 1:T) {
	      V[i, t] <- theta[1] + theta[2] * t/T + theta[3] * (t/T)^2
    }
}

# Generate the dependent variable Yit
beta <- c(0.5, 0.5)
Y <- array(0, dim = c(T, n))
for (i in 1:n) {
	  for (t in 1:T) {
		      Y[t, i] <- sum(beta * X[t, i, ]) + V[i, t] + rnorm(1)
  }
}

# Reshape Y for the KSS function
Y_mat <- matrix(Y, nrow = T, ncol = n)

# Apply the KSS function
results <- KSS(Y_mat ~ X[,,1] + X[,,2], additive.effects = "individual")
# summary(results)
# DGP5: Individual effects with large temporal variations

# Initialize the matrix for individual effects
V5 <- matrix(0, nrow = n, ncol = T)

# Generate time-varying individual effects vi(t) for DGP5
for (i in 1:n) {
	  upsilon1 <- rnorm(1)  # Coefficient for the first trigonometric function
  upsilon2 <- rnorm(1)  # Coefficient for the second trigonometric function
    for (t in 1:T) {
	        g1t <- sin(pi * t / 4)  # First trigonometric function
      g2t <- cos(pi * t / 4)  # Second trigonometric function
          V5[i, t] <- upsilon1 * g1t + upsilon2 * g2t  # Combine the effects
        }
}

# Generate the dependent variable Yit for DGP5
Y5 <- array(0, dim = c(T, n))
for (i in 1:n) {
	  for (t in 1:T) {
		      Y5[t, i] <- sum(beta * X[t, i, ]) + V5[i, t] + rnorm(1)
  }
}

# Reshape Y for the KSS function
Y5_mat <- matrix(Y5, nrow = T, ncol = n)

# Apply the KSS function
results5 <- KSS(Y5_mat ~ X[,,1] + X[,,2], additive.effects = "individual")
KSS.summary <- summary(results5)
print(KSS.summary)



