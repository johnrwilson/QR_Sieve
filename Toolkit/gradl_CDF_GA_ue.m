function [cdf_lf,cdf_g] = gradl_CDF_GA_ue(para_GA, taugrid, nmixtures, Y, X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradl_CDF_GA_ue
% Calculates the likelihood value and gradient for an uneven grid
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get number of covariates
    [nocovar,~] = size(X);

    % Part A) Convert para_GA to para_Pre
    % A-1) Convert step increase in para_GA to slope in para_Pre used in gradl_CDF_ue 
    step_size_ue = taugrid([2:end]) - taugrid([1:end-1]);
    ncbeta = length(taugrid);
    para = vec2mat(para_GA(1:nocovar*ncbeta), ncbeta);
    para(:,[2:end]) = para(:,[2:end])./(ones(nocovar,1)*step_size_ue);
    
    % A-2) Change the order of parameters from the GA setup to the new setup
    % Step 1: set the distributional parameters
    para_Pre = para_GA;
    para_Pre([1:nocovar*ncbeta]) = nan;
        
    % Step 2: set the constants     
    para_Pre([1:nocovar]) = para(:,1)';
    
    % Step 3: set the slope parameters
    para_temp = para(:,[2:end])';
    para_Pre([(nocovar+1) : (nocovar*ncbeta)]) = para_temp(:)';


    % Part B) Evaluate the log likelihood value and gradients
    [cdf_lf,cdf_g_new] = gradl_CDF_ue(para_Pre, taugrid, nmixtures, Y, X);

    
    % Part C) Set the gradients
    % Step 1: set the gradients for the distributional parameters
    cdf_g = cdf_g_new;
    cdf_g([1:nocovar*ncbeta]) = nan;
    
    
    % Step 2: set the constant and slope parameters
    g_increment = vec2mat(cdf_g_new(nocovar+1:nocovar*ncbeta), ncbeta-1)./(ones(nocovar,1)*step_size_ue);
    g_c = cdf_g_new(1:nocovar)';
    g_matrix = [g_c, g_increment]';
    
    cdf_g([1:nocovar*ncbeta]) = g_matrix(:)';
    
