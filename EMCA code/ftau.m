% f(y-xb(t))
function value = ftau(r_matrix, ParameterDist, nmixtures);
    lambda = ParameterDist(1:nmixtures)';
    mu = ParameterDist(nmixtures+1: 2*nmixtures)';
    sigma = ParameterDist(2*nmixtures+1: 3*nmixtures)';
    
    z = 0*r_matrix;
    for j = 1:nmixtures
        z = z+lambda(j)/sigma(j)*normpdf((r_matrix-mu(j))/sigma(j));
    end
    value = z;
end

