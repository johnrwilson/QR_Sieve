% In this version DGP has x1 and x1^2, but in the regression, we only
% include x1
% Haoyang Liu
% 4/27/2019

clear;
close all;

delete(gcp('nocreate'))
%parpool(4)

% A) Definition of constants
% Number of iterations (or runs of simulations)
n_repeatsample = 40;
ntau = 15;
n_WLS_iter = 40;


% taugrid is the old tau grid for quantile regression
% taugrid_midpoint is the mid point of tau segments
% taugrid_ue is the new uneven grid
[taugrid, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau);

% Part B) Preallocation of result variables
nsample = 100000;
ncovar = 2;
%number of mixture components, using mixture of normals
nmixtures=3;
%number of true components
nmixtures_truth = 3;
% nvars is the total number of parameters
nvars = ncovar*ntau+3*nmixtures-2;

% Preallocation for both quantile regression and MLE
% There are 3 types of estimates: qreg, qreg start (piecewise
% constant MLE), pc start (piecewise linear MLE)
recorder_qreg_repeat = nan(n_repeatsample,ncovar*ntau);
recorder_qreg_repeat_reshape = nan(n_repeatsample,ncovar,ntau);

% This extra set of qreg estimates are for y^* instead of y. 
recorder_qreg_repeat_1 = nan(n_repeatsample,ncovar*ntau);
recorder_qreg_repeat_reshape_1 = nan(n_repeatsample,ncovar,ntau);

recorder_WLS_repeat = nan(n_repeatsample,ncovar*ntau);
recorder_WLS_repeat_sorted = nan(n_repeatsample,ncovar,ntau);

% piecewise constant result
recorder_WLS_start_repeat = nan(n_repeatsample,nvars);
fval_recorder_WLS_start_repeat = nan(1,n_repeatsample);
exit_recorder_WLS_start_repeat = nan(1,n_repeatsample);
recorder_WLS_start_repeat_sorted = nan(n_repeatsample,ncovar,ntau);

% Mean of covariates
X_mean_repeat = nan(n_repeatsample,ncovar);

%Define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower=[zeros(1,(ncovar*ntau))-1, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(ncovar*ntau))+3), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];

%  Constants
b=1;
A = zeros(1,nvars);
A(ncovar*ntau+1) = 1;
A(ncovar*ntau+2) = 1;

% Part C) true coefficients
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

%true parameters for the error term:
sv_true=[lambda1_true,lambda2_true,mu1_true,mu2_true,1,1,1];


for j = [1:n_repeatsample]
    rng(j)
    j
    % Option for fmincon
    options=optimoptions(@fmincon,'GradObj','on');
    opts = optimoptions('ga', 'MaxGenerations',500,'PopulationSize',500);
    
    % D) Simulate the data
    tau_simu = rand(1,nsample);
    beta0_simu = b0(tau_simu);
    beta1_simu = b1(tau_simu);
    beta2_simu = b2(tau_simu);
    
    x1r = 2*2.34*(rand(1,nsample));
    x2r = 2*2.34*(x1r/(2*2.34)).^2;

    
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

    X = [ones(1,nsample); x1r]';
    
    X_mean = mean(X);
    X_mean_repeat(j,:) = X_mean;

    % Part E) qreg and WLS
    [fit] = quantlsfVector(X,y,taugrid);
    fit_1 = reshape(fit, ncovar*ntau,1);
    recorder_qreg = fit_1';
    recorder_qreg_repeat(j,:) = recorder_qreg;
    recorder_qreg_repeat_reshape(j,:,:) = fit';
    
    
    [fit_2nd] = quantlsfVector(X,y_s',taugrid);
    fit_1_2nd = reshape(fit, ncovar*ntau,1);
    recorder_qreg_1 = fit_1_2nd';
    recorder_qreg_repeat_1(j,:) = recorder_qreg_1;
    recorder_qreg_repeat_reshape_1(j,:,:) = fit_2nd';
    
    
    [recorder_WLS_repeat(j,:), recorder_WLS_repeat_sorted(j,:,:)] =  WLS_step(fit,X,y,n_WLS_iter);

    % Part F) MLE
    start = [recorder_WLS_repeat(j,:), para_dist_default];
    [fit_hat,fval,exitflag] = fmincon(@(x)gradllfCovarparavector(x, ntau, nsample, nmixtures,1,y, X),start,A,b,[],[],lower,upper,[],options);
    recorder_WLS_start_repeat(j,:) = fit_hat;
    fval_recorder_WLS_start_repeat(j) = fval;
    exit_recorder_WLS_start_repeat(j) = exitflag;


    % G) Piecewise linear MLE
    % G-1) Sort piecewise constant result

    recorder_WLS_start_reshape = reshape(fit_hat(1:ntau*ncovar),ntau,ncovar);
    [beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,ncovar);
    beta_WLS_start_sorted = beta_WLS_start_sorted';
    recorder_WLS_start_repeat_sorted(j,:,:) = beta_WLS_start_sorted;

  
end

recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted([1:j-1],:,:);
recorder_qreg_repeat_reshape = recorder_qreg_repeat_reshape([1:j-1],:,:);

clear X Y

%% 5) Save the result
save Ye_GA_3mixS_WLS_ue_ystar
return;
