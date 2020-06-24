%% 6) new error function introduced

function value = error_v2(theta, tgrid, t, dtau, x, y, nsample)
    value = repmat(y,length(t),1);
    k = 1;
    while k <= size(x, 1)
        value = value-repmat(x(k,:),length(t),1).*repmat(beta(theta(k,:), tgrid, t, dtau),1,nsample);
        k = k+1;
    end
end

