## functions for generating data

topic_poisson_generator <- function(n,p,k, seed = 0){
	## generate data X of size (p,n) with k hidden topics
	set.seed(seed)
	A = matrix(rnorm(p*k, 0, 1), nrow = p)
	A = exp(A)
	W = matrix(rnorm(k*n, 0, 1), nrow = k)
	W = exp(W)
	Lam = A %*% W
	X = matrix(rpois(n = p*n,lam = Lam), nrow = p)
	return(list(X = X, A = A, W = W))
}