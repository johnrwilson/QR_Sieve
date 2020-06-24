% This function reconstructs beta from alpha
function [beta] = reconstruct_beta(alpha)

    ntau = size(alpha,2);

    beta = alpha;
    %construct beta
    for j = 2:ntau
       beta(:,j) = beta(:,j-1) + alpha(:,j);
    end

    
end

