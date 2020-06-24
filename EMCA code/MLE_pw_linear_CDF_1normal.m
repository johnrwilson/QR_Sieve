% This codes incorporates Lei's new gradient codes using the CDF approach
% Write this code for just one normal error and see if 
% Haoyang Liu
% 3/3/2018

clear all;
close all;

%% 0) Option for the start value
% If truth_start == 1, use the truth as the start value
% If truth_start == 0, use qreg result as the start value
% If truth_start == -1, use the truth for c and beta, and default para dist values as the start value
% If truth_start == -2, use the truth for para dist, and qreg for beta and c as the start value
% If truth_start == -3, use qreg as the start value, except w_2
% If truth_start == -4, use the truth for parameters except c2
% If truth_start == -5, use the truth for parameters except b2(3)
truth_start = 0;

%% 1) Definition of constants
% Number of grid intervals. The number of grid points is one more than ntau
ntau = 33;

% In MLE, taugrid always covers the entire [0,1] interval. Thus the last
% grid point (1) is always omitted
taugrid = (0:(ntau-1))/(ntau);
taugrid_qreg = (0:(ntau))/(ntau);
epsilon = 0.00001;
taugrid_qreg(1) = epsilon + taugrid_qreg(1);
taugrid_qreg(end) =  -epsilon + taugrid_qreg(end);

% Number of iterations (or runs of simulations)
iter = 10;
nsample = 100000;

%number of mixture components, using mixture of normals
nmixtures=1;
%number of true components
nmixtures_truth = 1;

% nvars is the total number of parameters
% 3 for cs (constants)
% 3*ntau for beta (there are ntau intervals, and thus ntau slopes to estimate
% 3*nmixtures-2 for the distributional parameters
nvars=3+3*ntau+3*nmixtures-2;

% Option for fmincon
options=optimoptions(@fmincon,'display','iter','GradObj','on','UseParallel',true);%'DerivativeCheck','on');

% Default parameter for the distributions of measurement errors. They are
% used together with the result of WLS as the starting point for MLE
para_dist_default = [1];

% Preallocation for both quantile regression and MLE
% For qreg, tau grid is from 0 to 1. Thus there are ntau+1 points in the grid
recorder_qreg = nan(iter,3*ntau+3);
% For qreg start, the number of parameters is calculated above
recorder_qreg_start = nan(iter,nvars);
% The following are some one-dimensional parameters
fval_recorder_qreg_start = nan(iter,1);
exit_recorder_qreg_start = nan(iter,1);
iteration_recorder_qreg_start = nan(iter,1);
funcCount_recorder_qreg_start = nan(iter,1);
firstorderopt_recorder_qreg_start = nan(iter,1);

x1matrix = nan(iter,nsample);
x2matrix = nan(iter,nsample);
ytmatrix = nan(iter,nsample);

%true parameters for the error term:
sv_true=[1];

%define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower=[(zeros(1,(3*ntau+3))-20),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(3*ntau+3))+20), ones(1,nmixtures)*10];

% True qreg coefficients (level, not slope)
recorder_truth = [b0(taugrid_qreg), b1(taugrid_qreg), b2(taugrid_qreg)];

for j_iter = [1:iter]
%% 2) Simulate the data
    % Betas are simulated in a continuous range. So no need to change that. 
	tau_simu = rand(1,nsample);         
    beta0_simu = b0(tau_simu);
	beta1_simu = b1(tau_simu);
	beta2_simu = b2(tau_simu);
    
    x1r = 5*rand(1,nsample);
    x2r = 5*rand(1,nsample);
	% x1r = exp(normrnd(0,1,1,nsample));
    % x2r = exp(normrnd(0,1,1,nsample));
    
    x1matrix(j_iter,:) = x1r;
    x2matrix(j_iter,:) = x2r;

    y_ntemp = normrnd(0,sv_true,1,nsample);
    j_sample_begin = 1;
    [temp , y_n_index] = sort(rand(1,nsample));
    
    y_n = y_ntemp(y_n_index);
    y_s = beta0_simu + beta1_simu.*x1r + beta2_simu.*x2r;
    y = y_n+y_s;
    ytmatrix(j_iter,:) = y;        
    y = y';
    
%% 3) QREG
    X = [ones(1,nsample); x1r; x2r]';
    Y = repmat(y, [1 ntau]);
    
    [fit] = quantlsfVector(X,y,taugrid_qreg);        
    fit_1 = [squeeze(fit(:,1));squeeze(fit(:,2));squeeze(fit(:,3))];
    recorder_qreg(j_iter,:)=fit_1';
    
    
%% 4) MLE
      %
     display('reached MLE')  
     
     % Switch between different start values
     switch truth_start
         case -5 
             [c_start, beta_start] = pc_pl(recorder_truth, 3, (ntau));
             beta_start(33*2+5) = 1.8;
             start = [c_start, beta_start, sv_true];     
         case -4 
             [c_start, beta_start] = pc_pl(recorder_truth, 3, (ntau));
             c_start(2) = 0.8;
             start = [c_start, beta_start, sv_true];       
         case -3
             [c_start, beta_start] = pc_pl(squeeze(recorder_qreg(j_iter,:)), 3, (ntau));
             para_dist_default(2) = 0.25;
             start = [c_start, beta_start, para_dist_default];  
         case -2 
             [c_start, beta_start] = pc_pl(squeeze(recorder_qreg(j_iter,:)), 3, (ntau));
             start = [c_start, beta_start, sv_true]; 
         case -1 
             [c_start, beta_start] = pc_pl(recorder_truth, 3, (ntau));
             start = [c_start, beta_start, para_dist_default];     
         case 1 
             [c_start, beta_start] = pc_pl(recorder_truth, 3, (ntau));
             start = [c_start, beta_start, sv_true];
         otherwise
             [c_start, beta_start] = pc_pl(squeeze(recorder_qreg(j_iter,:)), 3, (ntau));
             start = [c_start, beta_start, para_dist_default];  
     end
     
     [start] = changestart( start,lower,upper );
     llf_start = gradl_CDF_Lei(start,  taugrid, nmixtures, y', X');

     if llf_start == Inf
            exit_recorder_qreg_start(j_iter) = -10;
     else
       % exit_recorder_qreg_start(j_iter) = 10;
       display('reached MLE fmincon')
      [fit_hat,fval,exitflag,output] = fmincon(@(x)gradl_CDF_Lei(x,  taugrid, nmixtures, y', X'),start,[],[],[],[],lower,upper,[],options);
     
      recorder_qreg_start(j_iter,:) = fit_hat;
      fval_recorder_qreg_start(j_iter) = fval;
      exit_recorder_qreg_start(j_iter) = exitflag;
      iteration_recorder_qreg_start(j_iter) = output.iterations;
      funcCount_recorder_qreg_start(j_iter) = output.funcCount;
      firstorderopt_recorder_qreg_start(j_iter) = output.firstorderopt;
     end

    
end
clear X Y


%% 5) Save the result
switch truth_start
    case -5
        save (sprintf('pl_CDF_truth_b2_5_%d',ntau))
    case -4
        save (sprintf('pl_CDF_truth_c2_%d_1norm',ntau))
    case -3
        save (sprintf('pl_CDF_qreg_3_%d_1norm',ntau))
    case -2
        save (sprintf('pl_CDF_qreg_para_%d_1norm',ntau))
    case -1
        save (sprintf('pl_CDF_truth_para_%d_1norm',ntau))
    case 1
        save (sprintf('pl_CDF_truth_%d_1norm',ntau))
    otherwise
        save (sprintf('pl_CDF_qreg_%d_1norm',ntau))
end
     

