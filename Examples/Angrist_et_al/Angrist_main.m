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
% data = [data; csvread('Angrist00.csv',1,0)];
% data = [data; csvread('Angrist90.csv',1,0)];
% data = [data; csvread('Angrist80.csv',1,0)];

ntau = 15;
% Number of WLS iterations
n_WLS_iter = 40;
nmixtures = 3;

% Part B) Load original data & keep WHITE observations
[X, y] = white_subsample(data);
clear data

% Calculate the number of observations in the subsamples
[nsample_orig, ncovar] = size(X);

nvars = ncovar*ntau+3*nmixtures-2;

% Upper and lower bounds for parameters in the MLE estimation
% lower=[zeros(1,(ncovar*ntau))-10000, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
lower = [-10000, .001, -10, .01];
% upper=[(zeros(1,(ncovar*ntau))+10000), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
upper = [10000, 1, 10, 10];

% Constants for MLE estimation, forcing the sum of the first two weights
% for mixtures not exceeding 1
b=1;
A = zeros(1,nvars);
A(ncovar*ntau+1) = 1;
A(ncovar*ntau+2) = 1;

do_mle = true;

n_batches = 50;
n_epochs = 500;
learning_rate = 0.00001;
decay = .999;

optimizer_settings = {'SGD', n_batches, n_epochs, learning_rate, decay, true};

[betas_mle, fit_hat, betas_bootstrap, fit_bootstrap] = QR_sieve(X, y, ntau, n_WLS_iter, upper, lower, nmixtures, A, b, do_mle, optimizer_settings, 1);