
function value = fun_beta_1(theta, tgrid, t, dtau, x, x2, y, taulower, ParameterDist, nmixtures)
    erroreachsample = error_Lei(theta, tgrid, t, dtau, x, y);
    value = fsum(erroreachsample, ParameterDist, nmixtures).*x2.*(t-taulower);
end
