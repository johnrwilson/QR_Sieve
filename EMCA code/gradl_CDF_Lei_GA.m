% This function 
% 1) Its parameter definition is the same as gradl_CDF_GA
% In other words, para_GA is the same as para_GA in gradl_CDF_GA
% 2) The function is implemented using gradl_CDF_Lei
function [cdf_lf,cdf_g] = gradl_CDF_Lei_GA(para_GA, taugrid, nmixtures, Y, X)
     [nocovar,nsample] = size(X);

    % Part A) Convert para_GA to para_Pre
    % A-1) Convert step increase in para_GA to slope in para_Pre used in gradl_CDF_Lei  
    stepsize = taugrid(3) - taugrid(2);
    ncbeta = length(taugrid);
    para = vec2mat(para_GA(1:nocovar*ncbeta), ncbeta);
    para(:,[2:end]) = para(:,[2:end])/stepsize;
    
    % A-2) Change the order of parameters from the GA setup to the Lei setup
    % Step 1: set the distributional parameters
    para_Pre = para_GA;
    para_Pre([1:nocovar*ncbeta]) = nan;
        
    % Step 2: set the constants     
    para_Pre([1:nocovar]) = para(:,1)';
    
    % Step 3: set the slope parameters
    para_temp = para(:,[2:end])';
    para_Pre([(nocovar+1) : (nocovar*ncbeta)]) = para_temp(:)';


    % Part B) Evaluate the log likelihood value and gradients
    taugrid = taugrid([1:(end-1)]);
    [cdf_lf,cdf_g_Lei] = gradl_CDF_Lei(para_Pre, taugrid, nmixtures, Y, X);

    
    % Part C) Set the gradients
    % Step 1: set the gradients for the distributional parameters
    cdf_g = cdf_g_Lei;
    cdf_g([1:nocovar*ncbeta]) = nan;
    
    
    % Step 2: set the constant and slope parameters
    g_slope = vec2mat(cdf_g_Lei(nocovar+1:nocovar*ncbeta), ncbeta-1)/stepsize;   
    g_c = cdf_g_Lei(1:nocovar)';
    g_matrix = [g_c, g_slope]';
    
    cdf_g([1:nocovar*ncbeta]) = g_matrix(:)';
    
