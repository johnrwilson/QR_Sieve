function [fit_1_temp] = construct_pl_start(beta_WLS_start_sorted, ncovar, ntau)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct_pl_start
% Construct the start value in the piecewise-linear estimator
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


beta_WLS_start_sorted(:, 1) = beta_WLS_start_sorted(:, 1) - 1/2*(beta_WLS_start_sorted(:, 2) - beta_WLS_start_sorted(:, 1));
beta_WLS_start_sorted(:, end) = beta_WLS_start_sorted(:, end) + 1/2*(beta_WLS_start_sorted(:, end) - beta_WLS_start_sorted(:, end-1));


fit_1_temp = nan(1,ncovar*(ntau));
for j_covar = [1:ncovar]
    % The start and the end index
    para_start = 2 + (j_covar-1)*(ntau);
    para_end = j_covar*(ntau);
    % The first constant
    fit_1_temp(1, para_start-1) = beta_WLS_start_sorted(j_covar, 1);
    % All increments after the first constant.
    fit_1_temp(1, [para_start : para_end]) = beta_WLS_start_sorted(j_covar, [2:end]) - beta_WLS_start_sorted(j_covar, [1:end-1]);
end
