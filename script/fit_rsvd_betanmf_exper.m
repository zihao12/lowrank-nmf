% experiment on rsvd + betanmf_mu (k topics)
% read Xhat computed from rsvd (D = 200)
% with desired d, reconstruct Xhat from rsvd output
% apply betanmf_exper on Xhat while keeping track of cost in both X and Xhat
% initialized with NNDSVD(X, k)

% hyperparameters
%dataname = 'test';
d = 50;

datadir = '../bigdata';
readcountsfile   = join([dataname, '.csv']);
initfactorsfile  = join([dataname, '_factors_rough.csv']);
initloadingsfile = join([dataname, '_loadings_rough.csv']);

Ufile = join([dataname, '_u_rsvd.csv']);
Vfile = join([dataname, '_v_rsvd.csv']);
Dfile = join([dataname, '_d_rsvd.csv']);


% These variables specify the names of the output files.
outdir = '../bigdata';
factorsoutfile  = join([dataname, '_factors_rsvdbetanmf_d',num2str(d),'.csv']);
loadingsoutfile = join([dataname, '_loadings_rsvdbetanmf_d',num2str(d),'.csv']);
erroutfile = join([dataname, '_error_rsvdbetanmf_d',num2str(d),'.csv']);
% SET UP ENVIRONMENT
% ------------------
addpath ../code

% LOAD GTEX DATA
% --------------
fprintf('Loading GTEx data.\n');
readcountsfile = fullfile(datadir,readcountsfile);
counts = csvread(readcountsfile);
fprintf('Loaded %d x %d count matrix.\n',size(counts,1),size(counts,2));

% LOAD RSVD RESULT
% --------------
fprintf('Loading GTEx data.\n');
Ufile = fullfile(datadir,Ufile);
Vfile = fullfile(datadir,Vfile);
Dfile = fullfile(datadir,Dfile);
U = csvread(Ufile);
V = csvread(Vfile);
D = csvread(Dfile);
%% recontruct using first d components
countshat = U(:,1:d) * diag(D(1:d)) * V(:,1:d)';
%% project to nonnegative
countshat(countshat < 0) = 0;
fprintf('reconstruct %d x %d count matrix.\n',size(countshat,1),size(countshat,2));


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
[A B f1 f2] = betanmf_exper(counts,countshat,L0,F0',1e-06,1000);
timing = toc;
fprintf('Computation took %0.2f seconds.\n',timing);

%% err: cost on counts; cost on countshat
err = [f1(:), f2(:)];

% Convert the Poisson model parameters to the parameters for the
% multinomial model.
[F L] = poisson2multinom(B',A);

% Compute the multinomial likelihood for the nnmf solution.
f = loglikmultinom(counts,F,L);
fprintf('Multinomial likelihood on counts   : %0.12f\n',f);
f = loglikmultinom(countshat,F,L);
fprintf('Multinomial likelihood on countshat: %0.12f\n',f);

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
