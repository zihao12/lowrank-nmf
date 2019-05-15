%% function to compute loglikelihood

function [type,multinom_ll, pois_ll] = compute_loglik(X,A,W)
	%% for now we only have poisson mdoels
    [p,n] = size(X);
    type = "poisson";
    pois_ll = sum(sum(log(poisspdf(X,A * W))))/(n*p);
    Ahat = A * diag(sum(A,1).^-1);
    Wtilde =  diag(sum(A,1)) * W;
    What = Wtilde * diag(sum(Wtilde,1).^-1);
    theta = Ahat * What;
    multinom_ll = sum(sum(X .* log(theta + eps)));
end
% A = F
% W = L'
% 
% X = counts'
% type = "poisson"
% pois_ll = sum(sum(log(poisspdf(X,A * W))))
% Ahat = A * diag(sum(F,1).^-1) 
% Wtilde =  diag(sum(F,1)) * W
% What = Wtilde * diag(sum(Wtilde,1).^-1)
% 
% theta = Ahat * What
% multinom_ll = sum(sum(X .* log(theta + eps)))
