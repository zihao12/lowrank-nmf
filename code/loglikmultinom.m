% Compute the log-likelihood for the multinomial topic model. Input X is an
% n x p matrix of counts, F is a p x K matrix of "factors", and L is an n x
% K matrix of "loadings".
function y = loglikmultinom (X, F, L)
  y = sum(sum(X .* log(L * F' + eps)));
    
