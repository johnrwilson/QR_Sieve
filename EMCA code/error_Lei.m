
% erroreachsample
function value = error_Lei(theta, tgrid, t, dtau, x, y)
    value = y;
    k = 1;
    while k <= size(x, 1)
        value = value-x(k,:)*beta(theta(k,:), tgrid, t, dtau);
        k = k+1;
    end
end

