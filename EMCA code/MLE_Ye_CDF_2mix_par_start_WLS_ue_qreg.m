% QR only 
% Haoyang Liu
% 1/23/2019

clear;
close all;


% A) Definition of constants
% Number of iterations (or runs of simulations)
n_repeatsample = 100;
ntau = 99;
n_WLS_iter = 40;


% taugrid is the old tau grid for quantile regression
% taugrid_midpoint is the mid point of tau segments
% taugrid_ue is the new uneven grid
[taugrid, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau);

% Part B) Preallocation of result variables
nsample = 100000;
ncovar = 3;
%number of mixture components, using mixture of normals
nmixtures=3;

% nvars is the total number of parameters
nvars = ncovar*ntau+3*nmixtures-2;

% Preallocation for both quantile regression and MLE
% There are 3 types of estimates: qreg, qreg start (piecewise
% constant MLE), pc start (piecewise linear MLE)
recorder_qreg_repeat = nan(n_repeatsample,ncovar*ntau);
recorder_qreg_repeat_reshape = nan(n_repeatsample,ncovar,ntau);

% Mean of covariates
X_mean_repeat = nan(n_repeatsample,ncovar);

%number of true components
nmixtures_truth = 2;
% Mixing probability of each component
lambda1_true = 1/3;
lambda2_true = 2/3;
% Mean of each component
mu1_true = -4.36; 
mu2_true = 2.18;
% Preprocess the parameters
lambdavector_true = [lambda1_true, lambda2_true];
muvector_true = [mu1_true, mu2_true];

for j = [1:n_repeatsample]
    j
    % Option for fmincon
    options=optimoptions(@fmincon,'GradObj','on');
    opts = optimoptions('ga', 'MaxGenerations',500,'PopulationSize',500);
    
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
    [temp , y_n_index] = sort(rand(1,nsample));
    
    y_n = y_ntemp(y_n_index);
    y_s = beta0_simu + beta1_simu.*x1r + beta2_simu.*x2r;
    y = y_n+y_s;
    y = y';
    

    X = [ones(1,nsample); x1r; x2r]';
  
    
    X_mean = mean(X);
    X_mean_repeat(j,:) = X_mean;

    % Part E) qreg and WLS
    [fit] = quantlsfVector(X,y,taugrid);
    fit_1 = reshape(fit, ncovar*ntau,1);
    recorder_qreg = fit_1';
    recorder_qreg_repeat(j,:) = recorder_qreg;
    recorder_qreg_repeat_reshape(j,:,:) = fit';
  

end



clear X Y

%% 5) Save the result
save Ye_GA_2mix_3_WLS_ue_qreg
return;
