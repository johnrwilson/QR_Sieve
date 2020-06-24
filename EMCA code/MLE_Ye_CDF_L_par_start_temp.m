% This codes gives a start value to piecewise linear MLE simulations.
% This version uses GA, instead of fmincon
% Haoyang Liu
% 8/15/2018

clear all;
close all;



% 1) Definition of constants
% Number of grid intervals. The number of grid points is one more than ntau
ntau = 9;

% In MLE, taugrid always covers the entire [0,1] interval. Thus the last
% grid point (1) is always omitted
taugrid = (1:(ntau))/(ntau+1);

% Number of iterations (or runs of simulations)
iter = 16;
nsample = 100000;
ncovar = 3;

%number of mixture components, using mixture of normals
nmixtures=3;


% nvars is the total number of parameters
% 3 for cs (constants)
% 3*ntau for beta (there are ntau intervals, and thus ntau slopes to estimate
% 3*nmixtures-2 for the distributional parameters
nvars = ncovar*ntau+3*nmixtures-2;

% Preallocation for both quantile regression and MLE
% There are 3 types of estimates: qreg, qreg start (piecewise
% constant MLE), pc start (piecewise linear MLE)
recorder_qreg = nan(iter,3*ntau);
recorder_qreg_start = nan(iter,nvars);
recorder_pc_start = nan(iter,nvars);

x1matrix = nan(iter,nsample);
x2matrix = nan(iter,nsample);
ytmatrix = nan(iter,nsample);

% Constants
b=1;
A = zeros(1,nvars);
A(3*ntau+1) = 1;
A(3*ntau+2) = 1;

% true coefficients
beta0_true = b0(taugrid);
beta1_true = b1(taugrid);
beta2_true = b2(taugrid);


%Define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower=[zeros(1,(3*ntau))-1, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(3*ntau))+3), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
lower_pl=[zeros(1,(3*ntau))-0.1, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper_pl=[(zeros(1,(3*ntau))+5), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];

for j_iter = [1:iter]
    j_iter
    % Option for fmincon
    options=optimoptions(@fmincon,'GradObj','on');%,'DerivativeCheck','on'
    opts = optimoptions('ga', 'MaxGenerations',500,'PopulationSize',500);%,'display','iter'
    
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
    
    y_n = randlap([1 nsample],1/2.29);
    
    y_s = beta0_simu + beta1_simu.*x1r + beta2_simu.*x2r;
    y = y_n+y_s;
    ytmatrix(j_iter,:) = y;
    y = y';
    
    %% 3) QREG
    X = [ones(1,nsample); x1r; x2r]';
    [fit] = quantlsfVector(X,y,taugrid);
    fit_1 = [squeeze(fit(:,1)); squeeze(fit(:,2)); squeeze(fit(:,3))];
    recorder_qreg(j_iter,:)=fit_1';
    true_beta = [beta0_true; beta1_true; beta2_true]';
    
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
    
    
    %% 5) MLE - piecewise linear
    % Sort piecewise constant result
    display('reached MLE (piecewise linear)');
    
    X_mean = mean(X)
    recorder_WLS_start_reshape = reshape(fit_hat(1:(ntau)*ncovar),(ntau),ncovar);
    [beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,(ntau),ncovar);
    beta_WLS_start_sorted = beta_WLS_start_sorted';
    
    % Construct the start value
    fit_1_temp = nan(1,ncovar*(ntau));
    for j_covar = [1:ncovar]
        % The start and the end index
        para_start = 2 + (j_covar-1)*(ntau);
        para_end = j_covar*(ntau);
        % The first constant
        fit_1_temp(1, para_start-1) = beta_WLS_start_sorted(j_covar, 1);
        % All increments after the first constant.
        fit_1_temp(1, [para_start : para_end]) = beta_WLS_start_sorted(j_covar, [2:end]) - beta_WLS_start_sorted(j_covar, [1:end-1]);
    end
    
    start = [fit_1_temp, fit_hat([(end - 3*nmixtures + 3) : end])];
    
    start = max(start, lower_pl+10^(-5));
    start = min(start, upper_pl-10^(-5));
    
    opts.InitialPopulationMatrix = start;
     
    taugrid_temp = [0:(ntau-1)]/(ntau-1);
    if gradl_CDF_GA(start,  taugrid_temp, nmixtures, y', X') == Inf
        exit_recorder_pc_start(j_iter) = -10;
    else
        [fit_hat_pl,fval,exitflag,output] =  ga(@(x)gradl_CDF_GA(x,  taugrid_temp, nmixtures, y', X'), nvars, A, b,[],[],lower_pl,upper_pl,[],opts);
        recorder_pc_start(j_iter,:) = fit_hat_pl;
        fval_recorder_pc_start(j_iter) = fval;
        exit_recorder_pc_start(j_iter) = exitflag;
    end
    
    
end
clear X Y

%% 5) Save the result
save Ye_GA_L_3_par_start
return;
