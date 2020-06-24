
function value = fun_phi_u(theta, tgrid, t, dtau, x, y, mu, sigma)
    erroreachsample = error_Lei(theta, tgrid, t, dtau, x, y);
    value = normpdf((erroreachsample-mu)/sigma);
end

