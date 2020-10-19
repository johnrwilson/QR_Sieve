function [loss, gradient] = ols_grad(b, y_dat, X_dat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ols_grad
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    loss = (y_dat - X_dat*b)'*(y_dat - X_dat*b) / size(y_dat,1);
    gradient = -2 * X_dat' * (y_dat - X_dat*b) / size(y_dat,1);
end
