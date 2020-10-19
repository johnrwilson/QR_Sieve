function [lambda,mu,lambda3] = preprocesslambdamu(lambdapre,mupre)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocesslambdamu
% Preprocess the EIV parameters
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
    

