
function value = fun_u_phi_u(theta, tgrid, t, dtau, x, y, mu, sigma)
    erroreachsample = error_Lei(theta, tgrid, t, dtau, x, y);
    value = ((erroreachsample-mu)/sigma).*normpdf((erroreachsample-mu)/sigma);
end

