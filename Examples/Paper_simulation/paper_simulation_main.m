clear all;

addpath("../../Toolkit/")

tic;

% This set of code implements the paper's method, using the simulated
% example from section 4 of the paper. Note that for now, this code only
% runs serially (not in parallel) so the sample size is much smaller than
% for the results reported in the paper.

% TODO: Improve comments for clarity

rng(10, 'twister')

ntau = 15;
n_WLS_iter = 40;

% Part B) Preallocation of result variables
nsample = 100000;
% number of covariates for the regression
ncovar = 3;
%number of mixture components, using mixture of normals
nmixtures = 3;
%number of true components
nmixtures_truth = 3;
% nvars is the total number of parameters
nvars = ncovar*ntau+3*nmixtures-2;

% Part C) true coefficients
[taugrid, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau);
beta0_true = b0(taugrid);
beta1_true = b1(taugrid);
beta2_true = b2(taugrid);

% Mixing probability of each component
lambda1_true = 0.5;
lambda2_true = 0.25;

% Mean of each component
mu1_true = -3;
mu2_true = 2;

% Preprocess the parameters
[lambdavector_true,muvector_true,q_true] =  preprocesslambdamu([lambda1_true,lambda2_true],[mu1_true,mu2_true]);
lambda3_true = lambdavector_true(nmixtures_truth);
mu3_true = muvector_true(nmixtures_truth);

% D) Simulate the data
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

n_batches = 50;
n_epochs = 50;
learning_rate = 0.00001;
decay = .9999;
optimizer_settings = {'SGD', n_batches, n_epochs, learning_rate, decay, true};

% lower = [-1, .001, -10, .001];
% upper = [3, 1, 10, 10];

% [betas_mle, fit_hat, betas_bootstrap, fit_bootstrap] = QR_sieve(X, y, 4, [], [], [], lower, upper, ...
%     false, optimizer_settings, true);

[betas_mle, fit_hat] = QR_sieve(X, y, 1, [], [], 500, [], [], ...
    true, optimizer_settings, true);

%%      
% figure()
% plot(taugrid_ue, betas_mle(1,:))
% title("beta_1")
% 
% figure()
% plot(taugrid_ue, betas_mle(2,:))
% title("beta_2")
% 
% figure()
% plot(taugrid_ue, betas_mle(3,:))
% title("beta_3")
% betas_std = std(betas_bootstrap, 0, 3);
% 
% for i = 1:ncovar
%     figure()
%     hold on
%     plot(taugrid_ue, betas_mle(i,:), 'LineWidth', 2)
%     plot(taugrid_ue, betas_mle(i,:) + 1.96 * betas_std(i,:), '--', 'LineWidth', 2)
%     plot(taugrid_ue, betas_mle(i,:) - 1.96 * betas_std(i,:), '--', 'LineWidth', 2)
%     xlabel('$\tau$', 'interpreter', 'latex', 'FontSize', 16)
%     ylabel_text = sprintf("$\\beta_{%d}(\\tau)$", i);
%     ylabel(ylabel_text, 'interpreter', 'latex', 'FontSize', 16);
% end
% 
% lambdas_short = fit_hat(ncovar*ntau+1:ncovar*ntau+nmixtures-1);
% lambdas = [lambdas_short, 1 - sum(lambdas_short)];
% mus_short = fit_hat(ncovar*ntau+nmixtures:ncovar*ntau+2*nmixtures-2);
% mus = [mus_short, -sum(lambdas_short.*mus_short)/lambdas(end)];
% sigmas = fit_hat(end-2:end);
% 
% dom_min = min(mus - 3 * sigmas);
% dom_max = max(mus + 3 * sigmas);
% n_points = 100;
% dom = linspace(dom_min, dom_max, n_points);
% 
% density_y = zeros(1, n_points);
% 
% for i = 1:nmixtures
%     density_y = density_y + lambdas(i) * normpdf(dom, mus(i), sigmas(i));
% end
% 
% figure()
% plot(dom, density_y)
% xlabel('Measurement Error', 'FontSize', 16)
% ylabel('Density', 'FontSize', 16);