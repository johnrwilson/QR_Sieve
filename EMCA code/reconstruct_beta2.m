% This function reconstructs piecewise linear level coefficients from
% piecewise linear raw estimates
function [b1,b2,b3] = reconstruct_beta2(alpha)
   
    k = (length(alpha))/3;
    alpha1 = alpha(1:k);
    alpha2 = alpha((k+1):2*k);
    alpha3 = alpha((2*k+1):3*k);
    beta1 = zeros(1,k);
    beta2 = zeros(1,k);
    beta3 = zeros(1,k);
    beta1(1) = alpha1(1);
    beta2(1) = alpha2(1);
    beta3(1) = alpha3(1);
    
    %construct beta
    for i = 2:k
        beta1(i) = beta1(i-1)+alpha1(i);
        beta2(i) = beta2(i-1)+alpha2(i);
        beta3(i) = beta3(i-1)+alpha3(i);
    end
    
    b1 = beta1;
    b2 = beta2;
    b3 = beta3;
    
end

        

