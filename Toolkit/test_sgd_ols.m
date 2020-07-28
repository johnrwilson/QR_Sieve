clear all;

n_samples = 100000;
X = ones(n_samples, 3);
X(:,2) = rand(n_samples,1) * 10 - 5;
X(:,3) = rand(n_samples,1) * 3 + 5;

errors = normrnd(0, 2, n_samples,1);

beta_true = [4; .2; -1];

y = X*beta_true + errors;

OLS_results = (X'*X)\(X'*y);

start = [1;1;1];
% start = OLS_results;

n_batches = 100;
n_epochs = 2000;
learning_rate = .01;
decay = .9999;

ols_grad(OLS_results, y, X)

f = @(a,b,c) ols_grad(a,b,c);

verbose = true;

opt_sgd = sgd(f, start, y, X, n_batches, n_epochs, learning_rate, decay, verbose);

ols_grad(opt_sgd, y, X)
ols_grad(beta_true, y, X)
ols_grad(OLS_results, y, X)