

% function used to calcualte the integral of f(y-xb(t))
function value = fun(theta, tgrid, t, dtau, x, y, ParameterDist, nmixtures)
    erroreachsample = error_Lei(theta, tgrid, t, dtau, x, y);
    value = ftau(erroreachsample, ParameterDist, nmixtures);
end

