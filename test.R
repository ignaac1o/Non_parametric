
mu <- c(0, 0)
Sigma <- matrix(c(1, -0.71, -0.71, 2), nrow = 2, ncol = 2)
x=seq(-3,3,l=12)
xx=as.matrix(expand.grid(x,x))

proj_grad_norm <- function(x, mu, Sigma, s = 2) {
  
  # Gradient
  grad <- grad_norm(x = xx[1,], mu = mu, Sigma = Sigma)
  
  # Hessian
  Hess <- Hess_norm(x = xx[1,], mu = mu, Sigma = Sigma)
  
  # Eigenvectors Hessian
  eig_Hess <- t(apply(Hess, 3, function(A) {
    eigen(x = A, symmetric = TRUE)$vectors[, s]
  }))
  
  # Projected gradient
  proj_grad <- t(sapply(1:nrow(eig_Hess), function(i) {
    tcrossprod(eig_Hess[i, ]) %*% grad[i, ]
  }))
  
  # As an array
  return(proj_grad)
  
}




# Euler solution

mu <- c(0, 0)
Sigma <- matrix(c(1, -0.71, -0.71, 2), nrow = 2, ncol = 2)
ks::plotmixt(mus = mu, Sigmas = Sigma, props = 1, display = "filled.contour2",
             gridsize = rep(251, 2), xlim = c(-5, 5), ylim = c(-5, 5),
             cont = seq(0, 90, by = 10), col.fun = viridis::viridis)

x0 <- as.matrix(expand.grid(seq(-3, 3, l = 12), seq(-3, 3, l = 12)))
x <- matrix(NA, nrow = nrow(x0), ncol = 2)
N <- 500
h <- 0.5
phi <- matrix(nrow = N + 1, ncol = 2)
eps <- 1e-4
for (i in 1:nrow(x0)) {
  
  #Move along the flow curve
  phi[1, ] <- x0[i, ]
  for (t in 1:N) {
    
    # Euler update
    phi[t + 1, ] <- phi[t, ] + 
      h * proj_grad_norm(phi[t, ], mu = mu, Sigma = Sigma) /
      mvtnorm::dmvnorm(x = phi[t, ], mean = mu, sigma = Sigma)
    
    # Stopping criterion (to save computing time!)
    abs_tol <- max(abs(phi[t + 1, ] - phi[t, ]))
    rel_tol <- abs_tol / max(abs(phi[t, ]))
    if (abs_tol < eps | rel_tol < eps) break
    
  }
  
  # Save final point
  x[i, ] <- phi[t + 1, , drop = FALSE]
  
  # Plot lines and x0
  lines(phi[1:(t + 1), ], type = "l")
  points(x0[i, , drop = FALSE], pch = 19)
  
}














grad_norm <- function(x, mu, Sigma) {
  
  # Check dimensions
  x <- rbind(x)
  p <- length(mu)
  stopifnot(ncol(x) == p & nrow(Sigma) == p & ncol(Sigma) == p)
  
  # Gradient
  grad <- -mvtnorm::dmvnorm(x = x, mean = mu, sigma = Sigma) *
    t(t(x) - mu) %*% solve(Sigma)
  return(grad)
  
}


Hess_norm <- function(x, mu, Sigma) {
  
  # Check dimensions
  x <- rbind(x)
  p <- length(mu)
  stopifnot(ncol(x) == p & nrow(Sigma) == p & ncol(Sigma) == p)
  
  # Hessian
  Sigma_inv <- solve(Sigma)
  H <- apply(x, 1, function(y) {
    mvtnorm::dmvnorm(x = y, mean = mu, sigma = Sigma) *
      (Sigma_inv %*% tcrossprod(y - mu) %*% Sigma_inv - Sigma_inv)
  })
  
  # As an array
  return(array(data = c(H), dim = c(p, p, nrow(x))))
  
}


a=proj_grad_norm(xx[1,],)


m=seq(-3,3,l=20)
count=0
for (i in 1:length(m)) {
  if(m[i]<(-1)){count=count+1}
  t[i]=m[i+count]
  if(i+count==length(m))break
}




