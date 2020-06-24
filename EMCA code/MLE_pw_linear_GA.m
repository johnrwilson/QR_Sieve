% This codes incorporates Lei's new gradient codes using the CDF approach
% Haoyang Liu

% 1/20/2017

clear all;
close all;

%% 0) Option for the start value
% If truth_start == 1, use the truth as the start value
% If truth_start == 0, use qreg result as the start value
% If truth_start == -1, use the truth for c and beta, and default para dist values as the start value
% If truth_start == -2, use the truth for para dist, and qreg for beta and c as the start value
truth_start = 0;

%% 1) Definition of constants
% Number of grid intervals. The number of grid points is one more than ntau
ntau = 9;

% In MLE, taugrid always covers the entire [0,1] interval. Thus the last
% grid point (1) is always omitted
taugrid = (0:(ntau-1))/(ntau);
taugrid_qreg = (0:(ntau))/(ntau);
epsilon = 0.001;
taugrid_qreg(1) = epsilon + taugrid_qreg(1);
taugrid_qreg(end) =  -epsilon + taugrid_qreg(end);

% Number of iterations (or runs of simulations)
iter = 1;
nsample = 100000;

%number of mixture components, using mixture of normals
nmixtures=3;
%number of true components
nmixtures_truth = 3;

% nvars is the total number of parameters
% 3 for cs (constants)
% 3*ntau for beta (there are ntau intervals, and thus ntau slopes to estimate
% 3*nmixtures-2 for the distributional parameters
nvars=3+3*ntau+3*nmixtures-2;

% Option for fmincon
options=optimoptions(@ga,'display','iter');%'DerivativeCheck','on');

% Default parameter for the distributions of measurement errors. They are
% used together with the result of WLS as the starting point for MLE
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];

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

%define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower=[(zeros(1,(3*ntau+3))-20), (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(3*ntau+3))+20), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
lower_ParaDist = [(zeros(1,(nmixtures-1))+0.01),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
upper_ParaDist = [(zeros(1,(nmixtures-1))*0+0.97),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];

% True qreg coefficients (level, not slope)
recorder_truth = [b0(taugrid_qreg), b1(taugrid_qreg), b2(taugrid_qreg)];

for j_iter = [1:iter]
%% 2) Simulate the data
    % Betas are simulated in a continuous range. So no need to change that. 
	tau_simu = rand(1,nsample);         
    beta0_simu = b0(tau_simu);
	beta1_simu = b1(tau_simu);
	beta2_simu = b2(tau_simu);
    
	x1r = exp(normrnd(0,1,1,nsample));
    x2r = exp(normrnd(0,1,1,nsample));
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
    Y = repmat(y, [1 ntau]);
    
    [fit] = quantlsfVector(X,y,taugrid_qreg);        
    fit_1 = [squeeze(fit(:,1));squeeze(fit(:,2));squeeze(fit(:,3))];
    recorder_qreg(j_iter,:)=fit_1';
    

    recorder_qreg(j_iter,[1,1])=[1]';
    recorder_qreg(j_iter,[1,2:10])=[0]';
%% 4) MLE
      %
     display('reached MLE')  
     
     % Switch between different start values
     switch truth_start
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
     
     [ start ] = changestart( start,lower,upper );
     llf_start = gradl_CDF_Lei(start,  taugrid, nmixtures, y', X');

     if llf_start == Inf
            exit_recorder_qreg_start(j_iter) = -10;
     else
       % exit_recorder_qreg_start(j_iter) = 10;
       display('reached MLE ga')
      [fit_hat,fval,exitflag,output] = ga(@(x)gradl_CDF_Lei(x,  taugrid, nmixtures, y', X'),nvars,A,b,[],[],lower,upper,[],options);
     
      recorder_qreg_start(j_iter,:)=fit_hat;
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
    case -2
        save (sprintf('pl_GA_qreg_para_%d_%d',nsample,ntau))
    case -1
        save (sprintf('pl_GA_truth_para_%d_%d',nsample,ntau))
    case 1
        save (sprintf('pl_GA_truth_%d_%d',nsample,ntau))
    otherwise
        save (sprintf('pl_GA_qreg_%d_%d',nsample,ntau))
end
     

