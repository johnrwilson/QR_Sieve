% This code calcualte the gradients with discontinuous beta functions using
% Christopher's CDF method. I also compare it with the old integration
% method to ensure we get the right results.

% Lei Ma
% 1/26/2018


clear all

%% 1) Define the parameters

% tau grid
ntau = 9;
taugrid = (0:(ntau-1))/(ntau);
taustep = 1/ntau;

nmixtures = 3;
nocovar = 3;
nsample = 10000;

% Mixing probability of each component
lambda1_true = 1/2;
lambda2_true = 1/4;
lambdapre = [lambda1_true lambda2_true];
% Mean of each component
mu1_true = -3; 
mu2_true = 2;
mupre = [mu1_true mu2_true];
% sigma
sigmavector_true = [1 1 1];

% Preprocess the parameters
[lambdavector_true,muvector_true,lambda3] =  preprocesslambdamu([lambda1_true,lambda2_true],[mu1_true,mu2_true]);
lambda3_true = lambdavector_true(nmixtures);
mu3_true = muvector_true(nmixtures);
ParameterDist = [lambdavector_true, muvector_true, sigmavector_true];

% true beta coefficients
beta0_true = repmat(2, 1, length(taugrid));
beta1_true = repmat(1, 1, length(taugrid));
beta2_true = repmat(1, 1, length(taugrid));
% beta(tau) function parameters
para_beta = [beta0_true; beta1_true; beta2_true];

% true c 
c0_true = repmat(0,  1, length(taugrid));
c1_true = repmat(-3, 1, length(taugrid));
c2_true = repmat(2,  1, length(taugrid));
para_c = [c0_true; c1_true; c2_true];

%% 2) simulate the data
    j_iter = 1;
    tau_simu = rand(1,nsample);         
    beta0_simu = b0(tau_simu);
	beta1_simu = b1(tau_simu);
	beta2_simu = b2(tau_simu);
    
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

%% 3) Evaluate the gradients

% tau_0,..., tau_K
t = (0:ntau)/ntau;
erroreachsample = dcerror(para_beta, para_c, taugrid, t, X, Y, nsample);
xB_k = repmat(X(1,:),ntau,1).*repmat(beta0_true',1,nsample)+ ...
       repmat(X(2,:),ntau,1).*repmat(beta1_true',1,nsample)+ ...
       repmat(X(3,:),ntau,1).*repmat(beta2_true',1,nsample);

% F - Li, G - lambda, H - mu, P - sigma   
F = 0*erroreachsample;
G = zeros(length(t), nsample, nmixtures);
H = zeros(length(t), nsample, nmixtures);
P = zeros(length(t), nsample, nmixtures);
f = 0*erroreachsample;

for j = 1:nmixtures
    F = F+lambdavector_true(j)*normcdf(erroreachsample,muvector_true(j),sigmavector_true(j));
    G(:,:,j) = normcdf(erroreachsample,muvector_true(j),sigmavector_true(j));
    H(:,:,j) = -lambdavector_true(j)./sigmavector_true(j).*...
               normpdf((erroreachsample-muvector_true(j))./sigmavector_true(j),0,1);
    P(:,:,j) = -lambdavector_true(j)./(sigmavector_true(j).^2).*...
               (erroreachsample-muvector_true(j)).*normpdf(erroreachsample,muvector_true(j),sigmavector_true(j));
    f = f+lambdavector_true(j)./sigmavector_true(j)*normpdf(erroreachsample,muvector_true(j),sigmavector_true(j));
end

F_diff = zeros(ntau, nsample);
G_diff = zeros(length(t)-1, nsample, nmixtures);
H_diff = zeros(length(t)-1, nsample, nmixtures);
P_diff = zeros(length(t)-1, nsample, nmixtures);
f_diff = zeros(ntau, nsample);
C = 1./xB_k;
C_diff = zeros(ntau-1,nsample);

for k = 1:ntau
    F_diff(k,:) = F(k,:) - F(k+1,:);
    f_diff(k,:) = f(k,:) - f(k+1,:);
end

for k = 1:ntau-1
    C_diff(k,:) = C(k,:) - C(k+1,:);
end
C_diff = [C_diff; C(end,:)];

for j = 1:nmixtures
    for k = 1:ntau
        G_diff(k,:,j) = G(k,:,j)-G(k+1,:,j);
        H_diff(k,:,j) = H(k,:,j)-H(k+1,:,j);
        P_diff(k,:,j) = P(k,:,j)-P(k+1,:,j);
    end
end

% likelihood
F_value = F_diff./xB_k;
cdf_lf_eachsample = sum(F_value);
cdf_lf = sum(cdf_lf_eachsample)/nsample;

% lambda gradients
G_value = G_diff./repmat(xB_k,1,1,3);
cdf_glambdaeachsample = squeeze(sum(G_value, 1)./cdf_lf_eachsample);
cdf_glambdaall = sum(cdf_glambdaeachsample)/nsample;

% mu gradients
H_value = H_diff./repmat(xB_k,1,1,3);
cdf_gmueachsample = squeeze(sum(H_value, 1)./cdf_lf_eachsample);
cdf_gmuall = sum(cdf_gmueachsample)/nsample;

% sigma gradients
P_value = P_diff./repmat(xB_k,1,1,3);
cdf_gsigmaeachsample = squeeze(sum(P_value, 1)./cdf_lf_eachsample);
cdf_gsigmaall = sum(cdf_gsigmaeachsample)/nsample;

% beta gradients

B1_value = (f(2:end,:).*t(2:end)'.*xB_k-F_diff)./(xB_k.^2);
B2_value = f(2:end,:).*t(2:end)';
B2_value = B2_value(1:end-1,:);
B2_value = B2_value./xB_k(2:end,:);
B2_value = [B2_value; zeros(1,nsample)];

cdf_gbetaeachsample = repmat((B1_value-B2_value)./cdf_lf_eachsample, 1,1,nmixtures);
cdf_gbetaeachsample(:,:,1)= cdf_gbetaeachsample(:,:,1).*x0;
cdf_gbetaeachsample(:,:,2)= cdf_gbetaeachsample(:,:,2).*x1r;
cdf_gbetaeachsample(:,:,3)= cdf_gbetaeachsample(:,:,3).*x2r;
cdf_gbeta = -squeeze(sum(cdf_gbetaeachsample, 2))/nsample;

% c gradients
cdf_gceachsample = repmat(C_diff.*f(2:end,:)./cdf_lf_eachsample,1,1,nmixtures);
cdf_gceachsample(:,:,1)= cdf_gceachsample(:,:,1).*x0;
cdf_gceachsample(:,:,2)= cdf_gceachsample(:,:,2).*x1r;
cdf_gceachsample(:,:,3)= cdf_gceachsample(:,:,3).*x2r;
cdf_gc = -squeeze(sum(cdf_gceachsample, 2))/nsample;


%% 4) Evaluate the gradients using "integral"

int_lf_eachsample = integral(@(t)dcfun(para_beta, para_c, taugrid, t, X, Y, ParameterDist, nmixtures), 0, 1, 'ArrayValued', true);
int_lf = sum(int_lf_eachsample)/nsample;

phi_u = zeros(nsample, nmixtures);
u_phi_u = zeros(nsample, nmixtures);
u2_phi_u = zeros(nsample, nmixtures);

for j = 1:nmixtures
    phi_u(:, j)=integral(@(t)dcfun_phi_u(para_beta, para_c, taugrid, t, X, Y, muvector_true(j), sigmavector_true(j)), 0, 1,'ArrayValued', true);
    u_phi_u(:, j)=integral(@(t)dcfun_u_phi_u(para_beta, para_c, taugrid, t, X, Y, muvector_true(j), sigmavector_true(j)), 0, 1,'ArrayValued', true);
    u2_phi_u(:, j)=integral(@(t)dcfun_u2_phi_u(para_beta, para_c, taugrid, t, X, Y, muvector_true(j), sigmavector_true(j)), 0, 1,'ArrayValued', true);
end

glambdaeachsample = phi_u./repmat(sigmavector_true, nsample, 1)./repmat(int_lf_eachsample', 1, nmixtures);
gmueachsample = u_phi_u.*repmat(lambdavector_true, nsample,1)./repmat(sigmavector_true.^2, nsample,1)./repmat(int_lf_eachsample', 1, nmixtures);
gsigma1eachsample = -1*phi_u.*repmat(lambdavector_true, nsample,1)./repmat(sigmavector_true.^2, nsample,1)./repmat(int_lf_eachsample', 1, nmixtures);
gsigma2eachsample = u2_phi_u.*repmat(lambdavector_true, nsample,1)./repmat(sigmavector_true.^2, nsample,1)./repmat(int_lf_eachsample', 1, nmixtures);
gsigmaeachsample = gsigma1eachsample + gsigma2eachsample;

glambda = sum(glambdaeachsample)/nsample;
gmu = sum(gmueachsample)/nsample;
gsigma = sum(gsigmaeachsample)/nsample;

% test for beta gradients -- add a small perturbation to beta and c
epsilon = 0.001;
%beta0_true(2) = beta0_true(2) + epsilon;
para_beta2 = [beta0_true; beta1_true; beta2_true];
c2_true(4) = c2_true(4)+epsilon;
para_c2 = [c0_true; c1_true; c2_true];
t2 = (0:ntau)/ntau;
erroreachsample2 = dcerror(para_beta2, para_c2, taugrid, t, X, Y, nsample);
xB_k2 = repmat(X(1,:),ntau,1).*repmat(beta0_true',1,nsample)+ ...
        repmat(X(2,:),ntau,1).*repmat(beta1_true',1,nsample)+ ...
        repmat(X(3,:),ntau,1).*repmat(beta2_true',1,nsample);
F2 = 0*erroreachsample2;
for j = 1:nmixtures
    F2 = F2+lambdavector_true(j)*normcdf(erroreachsample2,muvector_true(j),sigmavector_true(j));
end

F2_diff = zeros(ntau, nsample);
for k = 1:ntau
    F2_diff(k,:) = F2(k,:) - F2(k+1,:);
end
% likelihood
F2_value = F2_diff./xB_k2;
cdf_lf_eachsample2 = sum(F2_value);
cdf_lf2 = sum(log(cdf_lf_eachsample2))/nsample;
cdf_lf = sum(log(cdf_lf_eachsample))/nsample;
g(x) = (cdf_lf-cdf_lf2)/epsilon;

%% A1) functions for CDF method

function value = dcbeta(para_beta,para_c,tgrid,t)
    t = t';
    tgrid = repmat(tgrid, length(t), 1);
	[tbelow,tbelow_index]=max((tgrid<t).*tgrid, [], 2);
    value = para_beta(tbelow_index)'.*t+para_c(tbelow_index)';
end

% error = y-x*beta(tau_k)
function value = dcerror(para_beta, para_c, tgrid, t, x, y, nsample)
    value = repmat(y,length(t),1);
    i = 1;
    while i <= size(x, 1)
        value = value-repmat(x(i,:),length(t),1).*repmat(dcbeta(para_beta(i,:),para_c(i,:),tgrid, t),1,nsample);
        i = i+1;
    end
end

function [result] = b0(r)
    result = 2.*r;
end

function [result] = b1(r)
	result = r-3;
end

function [result] = b2(r)
	result = r+2;
end

function [lambda,mu,lambda3] = preprocesslambdamu(lambdapre,mupre)

% Process lambda
lambda3 = 1-sum(lambdapre);
lambda = [lambdapre,lambda3];

% Process mu
if length(mupre)>=1
mu=[mupre,(0-(mupre*(lambdapre'))/(lambda3))];
return;
end
if isempty(mupre)
    mu = 0;
end
end

%% A2) functions for integral method

% erroreachsample
function value = dcerror2(para_beta, para_c, tgrid, t, x, y)
    value = y;
    k = 1;
    while k <= size(x, 1)
        value = value-x(k,:)*dcbeta(para_beta(k,:), para_c(k,:), tgrid, t);
        k = k+1;
    end
end

% function used to calcualte the integral of f(y-xb(t))
function value = dcfun(para_beta, para_c, tgrid, t, x, y, ParameterDist, nmixtures)
    errorsample = dcerror2(para_beta, para_c, tgrid, t, x, y);
    value = ftau(errorsample, ParameterDist, nmixtures);
end

% f(y-xb(t))
function value = ftau(r_matrix, ParameterDist, nmixtures);
    lambda = ParameterDist(1:nmixtures)';
    mu = ParameterDist(nmixtures+1: 2*nmixtures)';
    sigma = ParameterDist(2*nmixtures+1: 3*nmixtures)';
    
    z = 0*r_matrix;
    for j = 1:nmixtures
        z = z+lambda(j)/sigma(j)*normpdf((r_matrix-mu(j))/sigma(j));
    end
    value = z;
end


function value = dcfun_phi_u(para_beta, para_c, tgrid, t, x, y, mu, sigma)
    erroreachsample = dcerror2(para_beta, para_c, tgrid, t, x, y);
    value = normpdf((erroreachsample-mu)/sigma);
end

function value = dcfun_u_phi_u(para_beta, para_c, tgrid, t, x, y, mu, sigma)
    erroreachsample = dcerror2(para_beta, para_c, tgrid, t, x, y);
    value = ((erroreachsample-mu)/sigma).*normpdf((erroreachsample-mu)/sigma);
end

function value = dcfun_u2_phi_u(para_beta, para_c, tgrid, t, x, y, mu, sigma)
    erroreachsample = dcerror2(para_beta, para_c, tgrid, t, x, y);
    value = (((erroreachsample-mu)/sigma).^2).* normpdf((erroreachsample-mu)/sigma);
end


