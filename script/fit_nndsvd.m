% fit data using NNDSVD
% arguments are: dataname & K & flag

% SCRIPT SETTINGS
% ---------------
% These variables specify the names of the input files.
%% dataname = "...." is input in batch file

%dataname = 'test';
%K = 20;
flag = 1;

datadir = '../bigdata';
readcountsfile   = join([dataname, '.csv']);

% These variables specify the names of the output files.
outdir = '../bigdata';
factorsoutfile  = join([dataname, '_factors_nndsvd_K',num2str(K),'.csv']);
loadingsoutfile = join([dataname, '_loadings_nndsvd_K',num2str(K),'.csv']);

% SET UP ENVIRONMENT
% ------------------
addpath ../code

% LOAD GTEX DATA
% --------------
fprintf('Loading GTEx data.\n');
readcountsfile = fullfile(datadir,readcountsfile);
counts = csvread(readcountsfile);
fprintf('Loaded %d x %d count matrix.\n',size(counts,1),size(counts,2));


% RUN NMF OPTIMIZATION METHOD
% ---------------------------
fprintf('Fitting with NNDSVD.\n')
tic;
[A,B] = NNDSVD(counts,K,flag);
timing = toc;
fprintf('Computation took %0.2f seconds.\n',timing);



% WRITE NNMF RESULTS TO FILE
% --------------------------
fprintf('Writing results to file.\n');
factorsoutfile  = fullfile(outdir,factorsoutfile);
loadingsoutfile = fullfile(outdir,loadingsoutfile);

csvwrite(factorsoutfile,B');
csvwrite(loadingsoutfile,A);

% SESSION INFO
% ------------
ver
