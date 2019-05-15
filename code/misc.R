
cost <- function(X,AB, e){
	return(sum(AB - X * log(AB + e)))
}

# Read a numeric (floating-point) matrix from a simple CSV file in
# which the row and column names are not given.
read.csv.matrix <- function (file) {
  out <- read_csv(file,col_names = FALSE,progress = FALSE,
                  col_types = cols(.default = col_double()))
  out <- as.matrix(out)
  rownames(out) <- NULL
  colnames(out) <- NULL
  return(out)
}

# Scale each column of A so that the entries in each column sum to 1;
# i.e., colSums(scale.cols(A)) should return a vector of ones.
scale.cols <- function (A)
  apply(A,2,function (x) x/sum(x))

# Convert the parameters (factors & loadings) for the Poisson model to
# the factors and loadings for the multinomial model. The return value
# "s" gives the Poisson rates for generating the "document" sizes.
poisson2multinom <- function (F, L) {
  L <- t(t(L) * colSums(F))
  s <- rowSums(L)
  L <- L / s
  F <- scale.cols(F)
  return(list(F = F,L = L,s = s))
}

# Compute the log-likelihood for the multinomial topic model. Input X
# is an n x p matrix of counts, F is a p x K matrix of "factors", and
# L is an n x K matrix of "loadings".
loglik.multinom <- function (X, F, L, e = .Machine$double.eps)
  sum(X * log(tcrossprod(L,F) + e))
