% Piecewise linear MLE on white subsamples using 
% This is the final version using a grid like [0.1:0.1:0.9] for quantile
% regression, and connecting the mid point values for pieceswise constant to 
% construct the start values for piecewise linear MLE
% Updated 12/9/2018

clear;
% close all;

addpath("../../Toolkit/")

data = csvread('census2010.csv',1,0);

% Part B) Load original data & keep white observations (need large N for MLE, see NBER version for effects by race x gender)
[X, y] = white_subsample(data);
clear data

% Upper and lower bounds for parameters in the MLE estimation
lower = [-10000, .001, -10, .01];
upper = [10000, 1, 10, 10];

n_batches = 50;
n_epochs = 500;
learning_rate = 0.00001;
decay = .999;

optimizer_settings = {'SGD', n_batches, n_epochs, learning_rate, decay, true};

[betas_mle, fit_hat, betas_bootstrap, fit_bootstrap] = sieve_mle(X, y, 1, [], [], [], lower, upper, optimizer_settings, true);
