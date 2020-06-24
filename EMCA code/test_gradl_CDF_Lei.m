% This code uses the CDF method to calculate the new gradients

clear all

%% 1) Input parameters
ntau = 9;
taugrid = (0:(ntau-1))/(ntau);
taustep = 1/ntau;

nmixtures = 3;
nocovar = 3;
nsample = 10000;

lambda1_true = 1/2;
lambda2_true = 1/4;
mu1_true = -3; 
mu2_true = 2;
sigmavector_true = [1 1 1];
% Preprocess the parameters
[lambdavector_true,muvector_true,lambda3] =  preprocesslambdamu([lambda1_true,lambda2_true],[mu1_true,mu2_true]);
ParameterDist = [lambdavector_true, muvector_true, sigmavector_true];

% NOTE: need to change the values when changing the beta(tau) functions
beta0_true = repmat(2, 1, length(taugrid));
beta1_true = repmat(1, 1, length(taugrid));
beta2_true = repmat(1, 1, length(taugrid));
cvector_true = [0 -3 2];

% beta(tau) function parameters
para = [beta0_true; beta1_true; beta2_true];
para = [cvector_true' para];

% transform the parameter format to the format Haoyang used
c = para(:,1);
betav = para(:,[2:end])';
beta_v = betav(:);
cbeta = [c', beta_v'];
cbeta = cbeta(:)';

%% 2) Simulate data for X and Y
    j_iter = 1;
    tau_simu = rand(1,nsample);         
    beta0_simu = b0_2(tau_simu);
	beta1_simu = b1_2(tau_simu);
	beta2_simu = b2_3(tau_simu);
    
    x1r = lognrnd(0,1,1,nsample);
    x2r = lognrnd(0,1,1,nsample);
    x1matrix(j_iter,:) = x1r;
    x2matrix(j_iter,:) = x2r;

    y_ntemp = zeros(1,nsample);
    j_sample_begin = 1;
    for j_mixture = [1 : nmixtures]
        if j_mixture == nmixtures
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
    
    x0 = ones(1, nsample);
    X = [x0; x1r; x2r];
    Y = y;


%% 3) calculate the gradients using the new CDF method

[cdf_lf, cdf_g] = gradl_CDF_Lei([cbeta, ParameterDist([1,2,4,5,7:9])], taugrid, nmixtures, Y, X);
% old method
[lf, g] = gradl_Lei([cbeta, ParameterDist([1,2,4,5,7:9])], taugrid, nmixtures, Y, X);




