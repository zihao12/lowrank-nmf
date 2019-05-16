% Compute the value of the objective function (the negative
% log-likelihood under the Poisson model) up to a normalizing
% constant.
function f = cost (X, AB, e)
	if nargin < 3
		e = 0;
	end	
  f = sum(sum(AB - X.*log(AB + e)));
