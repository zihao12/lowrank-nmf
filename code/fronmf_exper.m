% This function decomposes input matrix by nonnegative matrix factorization
% (NMF) based on beta-divergence criterion and multiplicative update rules.
%
% The input nonnegative data matrix X will be decomposed as X = A*B
% approximately, where matrices A and B are often called "basis matrix" and
% lee's mu update for NMF wiht Frobenius norm
%
function [A, B, f1, f2] = fronmf (X,Xhat, A, B, tol, maxiter, verbose)

  % Set zeros to small positive numbers to prevent numerical issues in
  % the updates below.
  [n p] = size(X);
  X     = max(X,eps);
  E     = ones(n,p);

  % Handle optional arguments.
  if nargin < 5
    tol = 1e-6;
  end
  if nargin < 6
    maxiter = 5000;
  end 
  if nargin < 7
    verbose = true;
  end

  % Compute the value of the objective function at the initial estimate
  % of the solution.
  f1    = zeros(maxiter + 1,1);
  f1(1) = costF(X,A*B);
  f2    = zeros(maxiter + 1,1);
  f2(1) = cost(Xhat,A*B);
  if verbose
    fprintf('iter obj (cost fn)  	obj_Xhat (cost fn)     max.diff\n');
    fprintf('---- ------------------- ------------------- --------\n');
  end

  % Repeat until the maximum number of iteration is reached, or until the
  % convergence criterion is met.
  for i = 1:maxiter
    A0     = A;
    B0     = B;
    [A B]  = update(Xhat,A,B);
    f1(i+1) = costF(X,A*B);
    f2(i+1) = costF(Xhat,A*B);
    d      = max(max(max(abs(A - A0))),...
                 max(max(abs(B - B0))));
    if verbose
      fprintf('%4d %+0.12e %+0.12e %0.2e\n',i,f1(i+1),f2(i+1),d);
    end
    if d < tol
      break
    end
  end
  f1 = f1(1:i+1);
  f2 = f2(1:i+1);

% This implements the multiplicative updates.
function [A, B] = update (X, A, B)
  A = A .* (X*B')./(A*(B*B'));
  A = max(A,eps);
  B = B .* (A'*X) ./ ((A'* A)*B);
  B = max(B,eps);
