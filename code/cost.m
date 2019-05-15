% Compute the value of the objective function (the negative
% log-likelihood under the Poisson model) up to a normalizing
% constant.
function f = cost (X, AB)
  f = sum(sum(AB - X.*log(AB + eps)));
