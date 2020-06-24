% This version
% This codes gives a start value to piecewise linear MLE simulations.
% This version uses GA, instead of fmincon
% Haoyang Liu
% 1/10/2019

clear;
close all;

%delete(gcp('nocreate'))
%parpool(25)


% A) Definition of constants
% Number of iterations (or runs of simulations)
n_repeatsample = 50;
ntau = 15;
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
recorder_WLS_repeat = nan(n_repeatsample,ncovar*ntau);

% piecewise constant result
recorder_WLS_start_repeat = nan(n_repeatsample,nvars);
fval_recorder_WLS_start_repeat = nan(1,n_repeatsample);
exit_recorder_WLS_start_repeat = nan(1,n_repeatsample);
recorder_WLS_start_repeat_sorted = nan(n_repeatsample,ncovar,ntau);



% Mean of covariates
X_mean_repeat = nan(n_repeatsample,ncovar);

%Define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower=[zeros(1,(3*ntau))-1, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(3*ntau))+3), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
lower_pl=[zeros(1,(3*ntau))-0.1, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper_pl=[(zeros(1,(3*ntau))+5), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];

% Constants
b=1;
A = zeros(1,nvars);
A(3*ntau+1) = 1;
A(3*ntau+2) = 1;

% Part C) true coefficients
beta0_true = b0(taugrid);
beta1_true = b1(taugrid);
beta2_true = b2(taugrid);


for j = [1:n_repeatsample]
    rng(j + 200)
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
    
    y_n = sqrt(3.5)*trnd(3,1,nsample);
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
    recorder_WLS_repeat(j,:) =  WLS_step(fit,X,y,n_WLS_iter);
    
    % Part F) MLE
    start = [recorder_WLS_repeat(j,:), para_dist_default];
    start = max(start,lower+10^(-5));
    start = min(start,upper-10^(-5));
    if gradllfCovarparavector(start, ntau, nsample, nmixtures,1, y, X) == Inf
        exit_recorder_WLS_start_repeat(j) = -10
    else
        [fit_hat,fval,exitflag] = fmincon(@(x)gradllfCovarparavector(x, ntau, nsample, nmixtures,1,y, X),start,A,b,[],[],lower,upper,[],options);
        recorder_WLS_start_repeat(j,:) = fit_hat;
        fval_recorder_WLS_start_repeat(j) = fval;
        exit_recorder_WLS_start_repeat(j) = exitflag;
    end
    
    j
    % G) Piecewise linear MLE
    % G-1) Sort piecewise constant result
    
    recorder_WLS_start_reshape = reshape(fit_hat(1:ntau*ncovar),ntau,ncovar);
    [beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,ncovar);
    beta_WLS_start_sorted = beta_WLS_start_sorted';
    recorder_WLS_start_repeat_sorted(j,:,:) = beta_WLS_start_sorted;
    
end



clear X Y

%% 5) Save the result
save Ye_GA_T_3_WLS_ue_5
return;
