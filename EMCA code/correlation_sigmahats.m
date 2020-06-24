% Load data version of the code
%
% Haoyang Liu
% 02/10/2015

clear all;
close all;

data=csvread('Angrist90.csv');
postfix = 'Angrist90';
% data = data([1:2000],:);
y=data(:,1);
X=[ones(size(y)), data(:,2:end)];

load WLS_60_Angrist90_99
[nsample, ncovar] = size(X);

% Use the results of 

%% 1) Definition of constants


% Preallocation
recorder_ParaDist = nan(1,(nvars-3*ntau));

% Constants
b=1;
A = zeros(1,nvars);
A(ncovar*ntau+1) = 1;
A(ncovar*ntau+2) = 1;


%% 2) WLS
betahat = fit';
Y = repmat(y, [1 ntau]);
ehat = Y-X*betahat;
sigmahat = sigmahat_WLS(60);
clear sigmahat_v
weight=ntau*normpdf(ehat./sigmahat)./repmat(sum(normpdf(ehat./sigmahat),2),[1 ntau]);
sigmahat = sigmahat_qreg;
weight_1=ntau*normpdf(ehat./sigmahat)./repmat(sum(normpdf(ehat./sigmahat),2),[1 ntau]);
    

weight_reshape = reshape(weight,86785*99,1);
weight_1_reshape = reshape(weight_1,86785*99,1);
corr(weight_reshape,weight_1_reshape)
    % j_WLS
