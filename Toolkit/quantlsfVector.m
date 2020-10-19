function [fit] = quantlsfVector(X,y,lp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantlsfVector
% Quantile regression on a vector of taus
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lv = length(lp);
[temp,lx] = size(X);
fit = NaN(lv,lx);
for jv = [1:lv] 
    [q  ] = lp(jv);
    [b] = quantlsf(X,y,q);
    fit(jv,:) = b;
end
