
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


# do log transform on count data X
# adopt the approach suggested by Abhishek Sarkar, and implemenetd in edgeR 
# X is (n_sample, n_feature) count matrix
# s = sum_j X_ij  (s_i is the sum for sample i)
# y_ij = ln( x_ij / (s_i/median(s)) + 1)  (edgeR uses mean instead of median)
# Y is our output
logtrans <- function(X){
  s = rowSums(X)
  w = s / median(s)
  Y = diag(1/w) %*% X + 1
  Y = log(Y)
  return(Y) 
}

# plot the histogram of residual and p value of fit using pearson test
# resid_ij = (x_ij - lam_ij)/(sqrt(lam_ij)) is normal under null, from which we compute p value
# two special cases:
# 1. x =0  , lam = 0. This will give NaN in R. I will set them to be 0
# 2. x ! = 0, lam = 0. This will give Inf in R, which is right

pval_resid_plot_pearson <- function(X, Lam,  resid.out.file, pval.out.file, main = " "){
  resids = (X-Lam)/sqrt(Lam)
  resids[is.na(resids)] = 0
  #resids[is.infinite(resids)] = sqrt(max(X))
  png(resid.out.file)
  hist(resids, breaks = 100, main = sprintf("Pearson Residual of %s", main))

  png(pval.out.file)
  pvals = 2*pnorm(-abs(resids))
  hist(pvals, breaks = 100, main = sprintf("Pearson pvalue of %s", main))
}

