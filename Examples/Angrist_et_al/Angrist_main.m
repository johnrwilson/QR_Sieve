% Piecewise linear MLE on white subsamples using 
% This is the final version using a grid like [0.1:0.1:0.9] for quantile
% regression, and connecting the mid point values for pieceswise constant to 
% construct the start values for piecewise linear MLE
% Updated 12/9/2018

clear;
% close all;

addpath("../../Toolkit/")

% data could be Angrist80 to Angrist00, or Jacob10
data=csvread('Jacob10.csv',1,0);
% data = [data; csvread('Angrist00.csv',1,0)];
% data = [data; csvread('Angrist90.csv',1,0)];
% data = [data; csvread('Angrist80.csv',1,0)];

% Part B) Load original data & keep WHITE observations
[X, y] = white_subsample(data);
clear data

% Upper and lower bounds for parameters in the MLE estimation
lower = [-10000, .001, -10, .01];
upper = [10000, 1, 10, 10];

n_batches = 50;
n_epochs = 100;
learning_rate = 0.00001;
decay = .999;

optimizer_settings = {'SGD', n_batches, n_epochs, learning_rate, decay, true};

[betas_mle, fit_hat, betas_bootstrap, fit_bootstrap] = QR_sieve(X, y, 1, [], [], [], [], [], [], optimizer_settings, true);