function [beta] = reconstruct_beta(alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstruct_beta
% reconstructs beta from alpha
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ntau = size(alpha,2);

    beta = alpha;
    %construct beta
    for j = 2:ntau
       beta(:,j) = beta(:,j-1) + alpha(:,j);
    end

    
end

