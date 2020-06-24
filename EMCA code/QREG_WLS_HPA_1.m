% This code analyzes the forecasting data used in Andreas's paper
%
% Haoyang Liu
% 12/30/2018

clear;
close all;
data=csvread('treatment.csv',1,0);

data = data(:,[1:5]);

% n_repeatsample could be 1
n_repeatsample = 1;
ntau = 15;
n_subsample = 1; % The number of observations in each subsample is calculated in Part B)
% Number of WLS iterations
n_WLS_iter = 30;

% taugrid is the old tau grid for quantile regression
% taugrid_midpoint is the mid point of tau segments
% taugrid_ue is the new uneven grid
% taugrid = (1:ntau)/(ntau+1);

taugrid_midpoint = (1/2*([0:(ntau-1)]  + [1:(ntau)]))/ntau;
taugrid_ue = taugrid_midpoint;
taugrid_ue(1) = 0;
taugrid_ue(end) = 1;


% Part B) Load original data & keep WHITE observations


% Keep
y_orig_1 = data(:,1);
X_orig_1 = [ones(1,length(y_orig_1));data(:,[2:end])']';


y_orig = y_orig_1((y_orig_1~=0&X_orig_1(:,2)==1));
X_orig = X_orig_1((y_orig_1~=0&X_orig_1(:,2)==1),:);
X_orig = X_orig(:,[1,4]);


[nsample_orig, ncovar] = size(X_orig);

[fit] = quantlsfVector(X_orig,y_orig,taugrid_midpoint);

for j_covar = [1:ncovar]
    figure;
    plot(taugrid_midpoint,fit(:,j_covar));
    title(sprintf('QREG coefficient for covariate %d',j_covar));
end



[recorder_WLS,recorder_WLS_sort] = WLS_step(fit,X_orig,y_orig,n_WLS_iter);


for j_covar = [1:ncovar]
    figure;
    plot(taugrid_midpoint,recorder_WLS_sort(:,j_covar));
    title(sprintf('WLS coefficient for covariate %d',j_covar));
end

return;
