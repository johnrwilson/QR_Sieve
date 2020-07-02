clear all;

addpath("../../Toolkit/")

% This set of code implements the paper's method, using the simulated
% example from section 4 of the paper. Note that for now, this code only
% runs serially (not in parallel) so the sample size is much smaller than
% for the results reported in the paper.

% TODO: Improve comments for clarity

ntau = 15;
n_WLS_iter = 40;

% Part B) Preallocation of result variables
nsample = 100;
% number of covariates for the regression
ncovar = 3;
%number of mixture components, using mixture of normals
nmixtures = 3;
%number of true components
nmixtures_truth = 3;
% nvars is the total number of parameters
nvars = ncovar*ntau+3*nmixtures-2;

%Define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
% beta,lambda,mu,sigma
lower=[zeros(1,(3*ntau))-1, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(3*ntau))+3), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];

% Constants
b=1;
A = zeros(1,nvars);
A(3*ntau+1) = 1;
A(3*ntau+2) = 1;

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

% Tell it to skip the MLE calculations for the sake of speed while testing
% the results.

% TODO: Make sure the code can handle the MLE calculations by running on
% server
n_batches = 10;
n_epochs = 5;
learning_rate = 0.001;
    
do_mle = false;

% betas_WLS_only = QR_sieve_test_sgd(X, y, ntau, n_WLS_iter, upper, lower, para_dist_default, A, b, do_mle, n_batches, n_epochs, learning_rate);

do_mle = true;

% betas_sgd = QR_sieve_test_sgd(X, y, ntau, n_WLS_iter, upper, lower, para_dist_default, A, b, do_mle, n_batches, n_epochs, learning_rate);

[betas_mle, fit_hat] = QR_sieve(X, y, ntau, n_WLS_iter, upper, lower, para_dist_default, A, b, do_mle);



% TODO: Plot the results in a way comparable to figures 1 and 2 from the
% paper

% if pl == 3
%     x1vector = [0:0.2:2*2.34];
%     x2vector = 2*2.34*(x1vector/(2*2.34)).^2;
%   
%     for j_tau = [2 8 14]
%         figure; 
%         hold on; 
%         plot(x1vector, beta_true(1,j_tau) + beta_true(2,j_tau)*x1vector + beta_true(3,j_tau)*x2vector,'b');
%         plot(x1vector, beta_qreg_full(1,j_tau) + beta_qreg_full(2,j_tau)*x1vector,'m');
%         plot(x1vector, beta_qreg_full_1(1,j_tau) + beta_qreg_full_1(2,j_tau)*x1vector,'m--'); 
%         plot(x1vector, beta_WLS_start_sorted_full(1,j_tau) + beta_WLS_start_sorted_full(2,j_tau)*x1vector,'k'); 
%         
%        lgd = legend('truth','quantile regression','quantile regression ystar','piecewise linear(fmincon)');  
%        c = lgd.Location;
%        lgd.Location = 'northwest';
%        title(sprintf('Q-tau(Y|X) tau = %0.1g',taugrid_ue(j_tau)))
%        print('-dpng','-r200',sprintf('qyx%d',j_tau));
%  
%     end
%  
% end

% plot(betas_mle(2,:))