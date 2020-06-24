% This codes tests the new gradient function gradl_CDF_Lei_GA 
% 1) If it matches gradl_CDF_GA (llf_Lei_GA == llf_GA)
% 2) If the analyatical gradient matches the empirical gradient (g_analytical == g_empirical)

% Haoyang Liu

% 3/14/2018

clear all;
close all;

%% 1) Definition of constants
% Number of grid intervals. The number of grid points is one more than ntau
ntau = 9;

% In MLE, taugrid always covers the entire [0,1] interval. Thus the last
% grid point (1) is always omitted
taugrid = (0:(ntau))/(ntau);
taugrid_qreg = (0:(ntau))/(ntau);
epsilon = 0.00001;
taugrid_qreg(1) = epsilon + taugrid_qreg(1);
taugrid_qreg(end) =  -epsilon + taugrid_qreg(end);

% Number of sample
nsample = 100;

%number of mixture components, using mixture of normals
nmixtures=3;
%number of true components
nmixtures_truth = 3;

% nvars is the total number of parameters
% 3 for cs (constants)
% 3*ntau for beta (there are ntau intervals, and thus ntau slopes to estimate
% 3*nmixtures-2 for the distributional parameters
nvars=3+3*ntau+3*nmixtures-2;

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

% True qreg coefficients (level, not slope)
recorder_truth = [b0(taugrid_qreg), b1(taugrid_qreg), b2(taugrid_qreg)];

%% 2) Simulate the data
    % Betas are simulated in a continuous range. So no need to change that. 
	tau_simu = rand(1,nsample);         
    beta0_simu = b0(tau_simu);
	beta1_simu = b1(tau_simu);
	beta2_simu = b2(tau_simu);
    
    x1r = 5*rand(1,nsample);
    x2r = 5*rand(1,nsample);
    X = [ones(1,nsample); x1r; x2r]';

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
    
    
%% 3) Create the truth 
     recorder_truth_temp = recorder_truth;
     recorder_truth_temp(1, [2:(ntau+1)]) = recorder_truth(1, [2:(ntau+1)]) - recorder_truth(1, [1:ntau]);
     recorder_truth_temp(1, [(ntau+3) : (2*ntau+2)]) = recorder_truth(1, [(ntau+3) : (2*ntau+2)]) - recorder_truth(1, [(ntau+2) : (2*ntau+1)]);
     recorder_truth_temp(1, [(2*ntau+4) : (3*ntau+3)]) = recorder_truth(1, [(2*ntau+4) : (3*ntau+3)])  - recorder_truth(1, [(2*ntau+3) : (3*ntau+2)]) ;
     
     
     start = [recorder_truth_temp, sv_true];
     
%% 4) Compare the empirical gradients and the analytical gradients      
     eps = 0.0001;
     
     g_analytical = zeros(1, nvars);
     g_empirical = zeros(1, nvars);
     g_empirical_GA = zeros(1, nvars);
     for j_var = [1:nvars]
         j_var
         start_temp = start;
         start_temp(j_var) = start_temp(j_var)*0.8;
         start_temp_eps = start_temp;
         start_temp_eps(j_var) = start_temp_eps(j_var) + eps;
         [llf, g]= gradl_CDF_Lei_GA(start_temp,  taugrid, nmixtures, y', X');
         [llf_1]= gradl_CDF_GA(start_temp,  taugrid, nmixtures, y', X');
         
         
         [llf_eps, g_1]= gradl_CDF_Lei_GA(start_temp_eps,  taugrid, nmixtures, y', X');
         
         
         [llf_1_eps]= gradl_CDF_GA(start_temp_eps,  taugrid, nmixtures, y', X');
         
         
         llf_Lei_GA(j_var) = llf;
         llf_GA(j_var) = llf_1;
         g_analytical(j_var) = g(j_var);
         g_empirical(j_var) = (llf_eps - llf)/(eps);
         g_empirical_GA(j_var) = (llf_1_eps - llf_1)/(eps);
         g_compare = [g_analytical;g_empirical;g_empirical_GA];
     end
     
     g_compare(4,:) = g_compare(1,:)./ g_compare(2,:);
     g_compare(5,:) = g_compare(1,:)./ g_compare(3,:);
     
     
     save test_gradl_CDF_Lei_GA_2
