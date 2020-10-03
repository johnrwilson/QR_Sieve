clear all;

addpath("../../Toolkit/")

tic;

% This set of code implements the paper's method, using the simulated
% example from section 4 of the paper. Note that the sample size used here 
% for dicactic purposes is is much smaller than for the results reported 
% in the paper, so the results will likely not be as smooth.

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate simulated data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed for testing different methods. This should be removed from the 
% final product.
rng(10, 'twister')

% Define some parameters for the model to run
ntau = 15;

nsample = 100000;
%number of true components
nmixtures_truth = 3;

% Mixing probability of each component
lambda1_true = 0.5;
lambda2_true = 0.25;

% Mean of each component
mu1_true = -3;
mu2_true = 2;

% True coefficients
[taugrid, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau);
beta0_true = b0(taugrid);
beta1_true = b1(taugrid);
beta2_true = b2(taugrid);

% Preprocess the parameters
[lambdavector_true,muvector_true,q_true] =  preprocesslambdamu([lambda1_true,lambda2_true],[mu1_true,mu2_true]);
lambda3_true = lambdavector_true(nmixtures_truth);
mu3_true = muvector_true(nmixtures_truth);

% Simulate the data
tau_simu = rand(1,nsample);
beta0_simu = b0(tau_simu);
beta1_simu = b1(tau_simu);
beta2_simu = b2(tau_simu);

x1r = exp(randn(1,nsample));
x2r = exp(randn(1,nsample));

y_ntemp = zeros(1,nsample);

j_sample_begin = 1;
for j_mixture = [1 : nmixtures_truth]
    if j_mixture == nmixtures_truth
        j_sample_end = nsample;
    else
        j_sample_end = min( ceil(sum(lambdavector_true(1:j_mixture))*nsample) , nsample);
    end
    y_ntemp([ j_sample_begin:j_sample_end]) = normrnd(muvector_true(j_mixture),1,1,[(j_sample_end - j_sample_begin +1)]);
    j_sample_begin = j_sample_end+1;
end
[~, y_n_index] = sort(rand(1,nsample));

y_n = y_ntemp(y_n_index);
y_s = beta0_simu + beta1_simu.*x1r + beta2_simu.*x2r;
y = y_n+y_s;
y = y';

X = [ones(1,nsample); x1r; x2r]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up for main toolkit %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Establish lower and upper bounds for each variable
% The order of the variables is [beta, lambda, mu, sigma].
lower = [-1, .001, -10, .001];
upper = [3, 1, 10, 10];

[betas_mle, fit_hat, betas_bootstrap, fit_bootstrap] = sieve_mle(X, y, 1, ...
                    [], [], [], lower, upper, [], true);