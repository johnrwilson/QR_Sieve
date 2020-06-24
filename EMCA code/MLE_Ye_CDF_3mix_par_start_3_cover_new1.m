% This code calculates the coverage probability
% This version uses two levels of loop to calculate the coverage
% probabilities
% Haoyang Liu
% 11/14/2018

clear;
close all;

delete(gcp('nocreate'))
parpool(25)

% 1) Definition of constants
% Number of grid intervals. The number of grid points is one more than ntau
ntau = 15;

% In MLE, taugrid always covers the entire [0,1] interval. Thus the last
% grid point (1) is always omitted
taugrid = (1:(ntau))/(ntau+1);

% Number of simulations
iter = 50;
% Number of bootstraps in each simulation
iter_each = 100;
nsample = 100000;
ncovar = 3;

%number of mixture components, using mixture of normals
nmixtures = 3;
%number of true components
nmixtures_truth = 3;

% nvars is the total number of parameters
% 3 for cs (constants)
% 3*ntau for beta (there are ntau intervals, and thus ntau slopes to estimate
% 3*nmixtures-2 for the distributional parameters
nvars = ncovar*ntau + 3*nmixtures - 2;

% Preallocation for both quantile regression and MLE
% There are 3 types of estimates: qreg, qreg start (piecewise
% constant MLE), pc start (piecewise linear MLE)

% The following nine matrices are for the first level of loop
recorder_qreg = nan(iter, ncovar*ntau);
recorder_qreg_start = nan(iter, nvars);
recorder_qreg_start_sorted = nan(iter, ncovar*ntau);
x1matrix = nan(iter, nsample);
x2matrix = nan(iter, nsample);
ytmatrix = nan(iter, nsample);
recorder_lambda_sorted = nan(iter, nmixtures);
recorder_mu_sorted = nan(iter, nmixtures);
recorder_sigma_sorted = nan(iter, nmixtures);

% The following three matrices are for the second level of loop
recorder_qreg_temp = nan(iter_each, ncovar*ntau);
recorder_qreg_start_temp = nan(iter_each, nvars);
recorder_qreg_start_sorted_temp = nan(iter_each, ncovar*ntau);
x1matrix_temp = nan(iter_each, nsample);
x2matrix_temp = nan(iter_each, nsample);
ytmatrix_temp = nan(iter_each, nsample);
mean_temp = nan(iter_each, ncovar);
recorder_lambda_sorted_temp = nan(iter_each, nmixtures);
recorder_mu_sorted_temp = nan(iter_each, nmixtures);
recorder_sigma_sorted_temp = nan(iter_each, nmixtures);

% The following three matrices are for the combined results from the two
% levels of loop
recorder_qreg_bootstrap = nan(iter, iter_each, ncovar*ntau);
recorder_qreg_start_bootstrap = nan(iter, iter_each, nvars);
recorder_qreg_start_bootstrap_sorted = nan(iter, iter_each, ncovar*ntau);
mean_bootstrap = nan(iter, iter_each, ncovar);
recorder_lambda_sorted_bootstrap = nan(iter, iter_each, nmixtures);
recorder_mu_sorted_bootstrap = nan(iter, iter_each, nmixtures);
recorder_sigma_sorted_bootstrap = nan(iter, iter_each, nmixtures);

% Constants
b = 1;
A = zeros(1,nvars);
A(ncovar*ntau+1) = 1;
A(ncovar*ntau+2) = 1;

% true coefficients
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
[lambdavector_true, muvector_true, q_true] =  preprocesslambdamu([lambda1_true,lambda2_true],[mu1_true,mu2_true]);
lambda3_true = lambdavector_true(nmixtures_truth);
mu3_true = muvector_true(nmixtures_truth);

%true parameters for the error term:
sv_true=[lambda1_true,lambda2_true,mu1_true,mu2_true,1,1,1];

%Define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower=[zeros(1,(3*ntau))-1, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(3*ntau))+3), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];

for j_iter = [1:iter]
    j_iter
    % Option for fmincon
    options = optimoptions(@fmincon,'GradObj','on');
    
    %% 2) Simulate the data
    % Betas are simulated in a continuous range. So no need to change that.
    tau_simu = rand(1,nsample);
    beta0_simu = b0(tau_simu);
    beta1_simu = b1(tau_simu);
    beta2_simu = b2(tau_simu);
    
    x1r = exp(randn(1,nsample));
    x2r = exp(randn(1,nsample));
    
    x1matrix(j_iter,:) = x1r;
    x2matrix(j_iter,:) = x2r;
    
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
    ytmatrix(j_iter,:) = y;
    y = y';
    
    %% 3) QREG
    X = [ones(1,nsample); x1r; x2r]';
    [fit] = quantlsfVector(X,y,taugrid);
    fit_1 = [squeeze(fit(:,1)); squeeze(fit(:,2)); squeeze(fit(:,3))];
    recorder_qreg(j_iter,:)=fit_1';
    
    
    %% 4) MLE - piecewise constant
    display('reached MLE (piecewise constant)');
    recorder_ParaDist = [[1/3,1/3],[-1,0],[1,1,1]];
    start = [fit_1',recorder_ParaDist];
    start = max(start,lower+10^(-5));
    start = min(start,upper-10^(-5));
    if gradllfCovarparavector(start, ntau, nsample, nmixtures,1, y, X) == Inf
        exit_recorder_qreg_start(j_iter) = -10;
    else
        [fit_hat,fval,exitflag] = fmincon(@(x)gradllfCovarparavector(x, ntau, nsample, nmixtures,1,y, X),start,A,b,[],[],lower,upper,[],options);
        fval_recorder_qreg_start(j_iter) = fval;
        recorder_qreg_start(j_iter,:) = fit_hat;
        exit_recorder_qreg_start(j_iter) = exitflag;
    end
    
    X_mean = [1,mean(x1r),mean(x2r)];
    recorder_qreg_start_sorted(j_iter,:) = sort_beta(X_mean,fit_hat,ncovar,ntau);
    
    
    [lambdasorted, musorted, sigmasorted] = sort_dist(fit_hat,ntau);
    recorder_lambda_sorted(j_iter,:) = lambdasorted;
    recorder_mu_sorted(j_iter,:) = musorted;
    recorder_sigma_sorted(j_iter,:) = sigmasorted;
    
    
    %% 5) Bootstrap
   parfor j_iter_each = [1:iter_each]
        
        
        %% 5-A) Simulate the data
        v_sample= randsample(nsample,nsample,true);
        
        x1r_temp = x1r(v_sample);
        x2r_temp = x2r(v_sample);
        y_temp = y(v_sample);
        
        
        %% 5-B) QREG
        X_temp = [ones(1,nsample); x1r_temp; x2r_temp]';
        [fit_temp] = quantlsfVector(X_temp,y_temp,taugrid);
        fit_1_temp = [squeeze(fit_temp(:,1)); squeeze(fit_temp(:,2)); squeeze(fit_temp(:,3))];
        recorder_qreg_temp(j_iter_each,:)=fit_1_temp';
        
        
        %% 5-C) MLE - piecewise constant
        start_temp = [fit_1_temp',recorder_ParaDist];
        start_temp = max(start_temp,lower+10^(-5));
        start_temp = min(start_temp,upper-10^(-5));
        if gradllfCovarparavector(start_temp, ntau, nsample, nmixtures,1, y_temp, X_temp) == Inf
            exit_recorder_qreg_start_temp(j_iter_each) = -10;
        else
            [fit_hat_temp,fval_temp,exitflag_temp] = fmincon(@(x)gradllfCovarparavector(x, ntau, nsample, nmixtures,1,y_temp, X_temp),start_temp,A,b,[],[],lower,upper,[],options);
            fval_recorder_qreg_start_temp(j_iter_each) = fval_temp;
            recorder_qreg_start_temp(j_iter_each,:) = fit_hat_temp;
            exit_recorder_qreg_start_temp(j_iter_each) = exitflag_temp;
        end
        
        
        X_mean_temp = [1,mean(x1r_temp),mean(x2r_temp)];
        mean_temp(j_iter_each,:) = X_mean_temp;
        recorder_qreg_start_sorted_temp(j_iter_each,:) = sort_beta(X_mean_temp,fit_hat_temp,ncovar,ntau);
        
        
        [lambdasorted_temp, musorted_temp, sigmasorted_temp] = sort_dist(fit_hat_temp,ntau);
        recorder_lambda_sorted_temp(j_iter_each,:) = lambdasorted_temp;
        recorder_mu_sorted_temp(j_iter_each,:) = musorted_temp;
        recorder_sigma_sorted_temp(j_iter_each,:) = sigmasorted_temp;
        
    end

    % 
    recorder_qreg_bootstrap(j_iter,:,:) = recorder_qreg_temp;
    recorder_qreg_start_bootstrap(j_iter,:,:) = recorder_qreg_start_temp;
    recorder_qreg_start_bootstrap_sorted(j_iter,:,:) = recorder_qreg_start_sorted_temp;
    mean_bootstrap(j_iter,:,:) = mean_temp;
    recorder_lambda_sorted_bootstrap(j_iter,:,:) = recorder_lambda_sorted_temp;
    recorder_mu_sorted_bootstrap(j_iter,:,:) = recorder_mu_sorted_temp;
    recorder_sigma_sorted_bootstrap(j_iter,:,:) = recorder_sigma_sorted_temp;
    
    
    save Ye_GA_3mix_3_par_start_3_cover_new_6
    
end
clear X Y

%% 5) Save the result
save Ye_GA_3mix_3_par_start_3_cover_new_6

