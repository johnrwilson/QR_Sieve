% Piecewise linear MLE on white subsamples using 
% This is the final version using a grid like [0.1:0.1:0.9] for quantile
% regression, and connecting the mid point values for pieceswise constant to 
% construct the start values for piecewise linear MLE
% Updated 12/9/2018

clear;
close all;

addpath("../../Toolkit/")

% Part A) Set parameters
% data could be Angrist80 to Angrist00, or Jacob10
data=csvread('Jacob10.csv',1,0);
postfix = 'Jacob10';

% n_repeatsample could be 1
n_repeatsample = 1;
ntau = 15;
n_subsample = 1; % The number of observations in each subsample is calculated in Part B)
% When n_subsample == 21, it runs the bootstrap procedure with the number
% of observations same as in the full sample. 

% Number of WLS iterations
n_WLS_iter = 40;
nmixtures = 3;

% taugrid is the old tau grid for quantile regression
% taugrid_midpoint is the mid point of tau segments
% taugrid_ue is the new uneven grid
[taugrid, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau);

% Part B) Load original data & keep WHITE observations
[X_orig, y_orig] = white_subsample(data);
clear data
% Calculate the number of observations in the subsamples
[nsample_orig, ncovar] = size(X_orig);
nsample = floor(1/n_subsample*nsample_orig);
% If n_subsample is 21, we do random sampling with replacement using the original number of observations
if n_subsample == 21
    nsample = nsample_orig;
end

nvars = ncovar*ntau+3*nmixtures-2;

% Upper and lower bounds for parameters in the MLE estimation
lower=[zeros(1,(ncovar*ntau))-10000, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(ncovar*ntau))+10000), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];

% Constants for MLE estimation, forcing the sum of the first two weights
% for mixtures not exceeding 1
b=1;
A = zeros(1,nvars);
A(ncovar*ntau+1) = 1;
A(ncovar*ntau+2) = 1;
    
if n_subsample ~= 1
    v_sample= randsample(nsample_orig,nsample,true);
else
    v_sample = [1:nsample_orig];
end
[v_sample_sorted] = sort(v_sample);
y = y_orig(v_sample);
X = X_orig(v_sample,:);

do_mle = false;

betas_WLS_only = QR_sieve(X, y, ntau, n_WLS_iter, upper, lower, para_dist_default, A, b, do_mle);