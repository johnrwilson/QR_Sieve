
function value = fun_u2_phi_u(theta, tgrid, t, dtau, x, y, mu, sigma)
    erroreachsample = error_Lei(theta, tgrid, t, dtau, x, y);
    value = (((erroreachsample-mu)/sigma).^2).* normpdf((erroreachsample-mu)/sigma);
end
