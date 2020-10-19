function [beta_WLS_start_sorted,I_WLS_start] = sortbeta_1(X_mean,beta_WLS_start,ntau,nx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sortbeta_1
% Sort knots based on the predicted values
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fval_tau_WLS_start = X_mean*beta_WLS_start';
[V,I_WLS_start] = sort(fval_tau_WLS_start);
for j_x = [1:nx]
    beta_WLS_start_sorted(:,j_x) = beta_WLS_start(I_WLS_start,j_x);
end
