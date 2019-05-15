% Convert the parameters (factors & loadings) for the Poisson model to the
% factors and loadings for the multinomial model. The return value "s" gives
% the Poisson rates for generating the "document" sizes.
function [F, L, s] = poisson2multinom (F, L)
  L = L .* sum(F);
  s = sum(L,2);
  L = L ./ s;
  F = F ./ sum(F);
