% This codes incorporates Lei's new gradient codes
% Haoyang Liu
% 11/17/2017

clear all;
close all;


%% 1) Definition of constants
% Number of grid points
ntau = 9;
% taugrid is a grid on [0,1]
taugrid = (0:(ntau-1))/(ntau);
taugrid_qreg = (0:(ntau))/(ntau);
epsilon = 0.002;
taugrid_qreg(1) = epsilon + taugrid_qreg(1);
taugrid_qreg(end) =  - epsilon + taugrid_qreg(end);

% Number of iterations (or runs of simulations)
iter = 10;
nsample = 10000;

%number of mixture components, using mixture of normals
nmixtures=3;
nmixtures_truth = 3;

% nvars is the total number of parameters
% 3 for cs, 3*ntau for betas, 3*nmixtures-2 for the distributional
% parameters
nvars=3+3*ntau+3*nmixtures-2;

% Option for fmincon
options=optimoptions(@fmincon,'display','iter','GradObj','on','UseParallel',true,'DerivativeCheck','on');

% Default parameter for the distributions of measurement errors. They are
% used together with the result of WLS as the starting point for MLE
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];

% Preallocation for both quantile regression and MLE
recorder_qreg = nan(iter,3*ntau+3);

recorder_qreg_start = nan(iter,nvars);
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
[lambdavector_true,muvector_true,q_true] =  preprocesslambdamu([lambda1_true,lambda2_true],[mu1_true,mu2_true]);
lambda3_true = lambdavector_true(nmixtures_truth);
mu3_true = muvector_true(nmixtures_truth);

%true parameter sv:
sv_true=[lambda1_true,lambda2_true,mu1_true,mu2_true,1,1,1];

%define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower=[(zeros(1,(3*ntau+3))-5), (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(3*ntau+3))+10), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
lower_ParaDist = [(zeros(1,(nmixtures-1))+0.01),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
upper_ParaDist = [(zeros(1,(nmixtures-1))*0+0.97),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];



for j_iter = [1:iter]
    display('got here')
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
    

 
%% 4) MLE
     %
      display('got MLE')
     [c_start, beta_start] = pc_pl(squeeze(recorder_qreg(j_iter,:)), 3, (ntau+1));
     start = [c_start, beta_start, para_dist_default];  
     [ start ] = changestart( start,lower,upper );
     llf_start = gradl_Lei(start,  taugrid, nmixtures, y', X');

     if llf_start == Inf
            exit_recorder_qreg_start(j_iter) = -10;
     else
       
        display('got MLE fmincon')
        [fit_hat,fval,exitflag,output] = fmincon(@(x)gradl_Lei(x,  taugrid, nmixtures, y', X'),start,A,b,[],[],lower,upper,[],options);
     
        
        
        recorder_qreg_start(j_iter,:)=fit_hat;
        fval_recorder_qreg_start(j_iter) = fval;
        exit_recorder_qreg_start(j_iter) = exitflag;
        iteration_recorder_qreg_start(j_iter) = output.iterations;
        funcCount_recorder_qreg_start(j_iter) = output.funcCount;
        firstorderopt_recorder_qreg_start(j_iter) = output.firstorderopt;
     end

    
end
