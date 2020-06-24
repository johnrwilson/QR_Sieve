% This codes gives a start value to simulations so that their problem is a
% bit easier.
% Haoyang Liu
% 8/13/2018

clear all;
close all;

% 1) Definition of constants
options=optimoptions(@fmincon,'GradObj','on','MaxFunEvals',30000,'MaxIter',10000);
opts = optimoptions('ga', 'MaxGenerations',5000,'PopulationSize',500,'UseParallel',true,'display','iter');

% Number of grid intervals. The number of grid points is one more than ntau
ntau = 10;
n_WLS_iter = 100;
% In MLE, taugrid always covers the entire [0,1] interval. Thus the last
% grid point (1) is always omitted
taugrid = (0:(ntau))/(ntau);
taugrid_qreg = (0:(ntau))/(ntau);
epsilon = 0.01;
taugrid_qreg(1) = epsilon + taugrid_qreg(1);
taugrid_qreg(end) =  -epsilon + taugrid_qreg(end);

% Number of iterations (or runs of simulations)
iter = 5;
nsample = 100000;
ncovar = 3;

%number of mixture components, using mixture of normals
nmixtures=3;
%number of true components
nmixtures_truth = 3;

% nvars is the total number of parameters
% 3 for cs (constants)
% 3*ntau for beta (there are ntau intervals, and thus ntau slopes to estimate
% 3*nmixtures-2 for the distributional parameters
nvars=3+3*ntau+3*nmixtures-2;

% Preallocation for both quantile regression and MLE
% There are 4 or 5 types of estimates: qreg, WLS, WLS start (piecewise
% constant MLE), pc start (piecewise linear MLE), and maybe fmincon start
% (piecewise linear MLE)
% For qreg, tau grid is from 0 to 1. Thus there are ntau+1 points in the grid
recorder_qreg = nan(iter,3*ntau+3);
recorder_WLS = nan(iter,3*ntau+3);
recorder_WLS_start = nan(iter,nvars);
recorder_pc_start = nan(iter,nvars);
recorder_fmincon_start = nan(iter,nvars);

x1matrix = nan(iter,nsample);
x2matrix = nan(iter,nsample);
ytmatrix = nan(iter,nsample);

% Constants
b=1;
A = zeros(1,nvars);
A(3*ntau+3+1) = 1;
A(3*ntau+3+2) = 1;

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


%Define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower_pl = [zeros(1,(ncovar*(ntau+1)))-10, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
upper_pl = [(zeros(1,(ncovar*(ntau+1)))+10), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
lower_pc = [zeros(1,(ncovar*(ntau+1)))-10000, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
upper_pc = [(zeros(1,(ncovar*(ntau+1)))+10000), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];

ncovar = 3;

for j_iter = [1:iter]
    j_iter
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
    Y = repmat(y, [1 (ntau+1)]);
    
    [fit] = quantlsfVector(X,y,taugrid_qreg);
    fit_1 = [squeeze(fit(:,1));squeeze(fit(:,2));squeeze(fit(:,3))];
    recorder_qreg(j_iter,:)=fit_1';
    
  
    
    %% 4) MLE - piecewise constant
    display('reached MLE (piecewise constant)');
    recorder_ParaDist = [[1/3,1/3],[-1,0],[1,1,1]];  
    start = [fit_1',recorder_ParaDist];
     
    [fit_hat,fval,exitflag,output] = fmincon(@(x)gradllfCovarparavector(x, ntau+1, nsample, nmixtures,1,y, X),start,A,b,[],[],lower_pc,upper_pc,[],options);
    recorder_WLS_start(j_iter,:) = fit_hat;
    
    
    %% 6) MLE - piecewise linear
    % Sort piecewise constant result
    X_mean = mean(X);
    recorder_WLS_start_reshape = reshape(fit_hat(1:(ntau+1)*ncovar),(ntau+1),ncovar);
    [beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,(ntau+1),ncovar);
    beta_WLS_start_sorted = beta_WLS_start_sorted'
    
    % Construct the start value
    fit_1_temp = nan(1,ncovar*(ntau+1));
    for j_covar = [1:ncovar]
        % The start and the end index
        para_start = 2 + (j_covar-1)*(ntau+1);
        para_end = j_covar*(ntau+1);
        
        % The first constant
        fit_1_temp(1, para_start-1) =  beta_WLS_start_sorted(j_covar, 1);
        % All increments after the first constant. 
        fit_1_temp(1, [para_start : para_end]) = beta_WLS_start_sorted(j_covar, [2:end]) - beta_WLS_start_sorted(j_covar, [1:end-1]);
    end
    
    start = [fit_1_temp, fit_hat([(end - 3*nmixtures + 3) : end])];
    [fit_hat,fval,exitflag,output] = fmincon(@(x)gradl_CDF_Lei_GA(x,  taugrid, nmixtures, y', X'), start, A, b,[],[],lower_pl,upper_pl,[],options);
    recorder_pc_start(j_iter,:) = fit_hat;
    fval_recorder_pc_start(j_iter) = fval;
    exit_recorder_pc_start(j_iter) = exitflag;
    iteration_recorder_pc_start(j_iter) = output.iterations;
    funcCount_recorder_pc_start_repeat(j_iter) = output.funcCount;
    firstorderopt_recorder_pc_start(j_iter) = output.firstorderopt;
    
    
end
clear X Y

%% 5) Save the result
save Ye_GA_3mix_3_par_start
return;
