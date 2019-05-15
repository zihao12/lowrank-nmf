% Fit a non-megative matrix factorization to the GTEx data using the
% multiplication update rules; see betanmf.m for more information about the
% algorithm, and see fit_gtex_betanmf.sbatch for SLURM settings used on the
% RCC cluster.

% SCRIPT SETTINGS
% ---------------
% These variables specify the names of the input files.
%% dataname = "...." is input in bacth file

datadir = '../bigdata';
readcountsfile   = join([dataname, '.csv']);
initfactorsfile  = join([dataname, '_factors_rough.csv']);
initloadingsfile = join([dataname, '_loadings_rough.csv']);

% These variables specify the names of the output files.
outdir = '../bigdata';
factorsoutfile  = join([dataname, '_factors_betanmf.csv']);
loadingsoutfile = join([dataname, '_loadings_betanmf.csv']);
erroutfile = join([dataname, '_error_betanmf.csv']);

% SET UP ENVIRONMENT
% ------------------
addpath ../code

% LOAD GTEX DATA
% --------------
fprintf('Loading GTEx data.\n');
readcountsfile = fullfile(datadir,readcountsfile);
counts = csvread(readcountsfile);
fprintf('Loaded %d x %d count matrix.\n',size(counts,1),size(counts,2));

% LOAD INITIAL ESTIMATES
% ----------------------
fprintf('Loading initial estimates of factors and loadings.\n');
initfactorsfile  = fullfile(datadir,initfactorsfile);
initloadingsfile = fullfile(datadir,initloadingsfile);
F0               = csvread(initfactorsfile);
L0               = csvread(initloadingsfile);
fprintf('Loaded %d x %d factors matrix, ',size(F0,1),size(F0,2));
fprintf('and %d x %d loadings matrix.\n',size(L0,1),size(L0,2));

% RUN NMF OPTIMIZATION METHOD
% ---------------------------
fprintf('Fitting Poisson topic model using betanmf.\n')
tic;
[A B err] = betanmf(counts,L0,F0',1e-06,1000);
timing = toc;
fprintf('Computation took %0.2f seconds.\n',timing);

% Convert the Poisson model parameters to the parameters for the
% multinomial model.
[F L] = poisson2multinom(B',A);

% Compute the multinomial likelihood for the nnmf solution.
f = loglikmultinom(counts,F,L);
fprintf('Multinomial likelihood at nnmf solution: %0.12f\n',f);

% WRITE NNMF RESULTS TO FILE
% --------------------------
fprintf('Writing results to file.\n');
factorsoutfile  = fullfile(outdir,factorsoutfile);
loadingsoutfile = fullfile(outdir,loadingsoutfile);
erroutfile = fullfile(outdir,erroutfile);

csvwrite(factorsoutfile,F);
csvwrite(loadingsoutfile,L);
dlmwrite(erroutfile,err,'precision', '%0.12e');

% SESSION INFO
% ------------
ver
