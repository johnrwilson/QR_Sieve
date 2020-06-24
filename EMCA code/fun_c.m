

function value = fun_c(theta, tgrid, t, dtau, x, x2, y, ParameterDist, nmixtures)
    erroreachsample = error_Lei(theta, tgrid, t, dtau, x, y);
    value = fsum(erroreachsample, ParameterDist, nmixtures).*x2;
end


