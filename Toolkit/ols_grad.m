function [loss, gradient] = ols_grad(b, y_dat, X_dat)
    loss = (y_dat - X_dat*b)'*(y_dat - X_dat*b) / size(y_dat,1);
    gradient = -2 * X_dat' * (y_dat - X_dat*b) / size(y_dat,1);
end