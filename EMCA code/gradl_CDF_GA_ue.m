%% gradl for uneven tau grid
% This code is exactly the same as gradl_CDF_GA. What needs edits is code
% gradl_CDF_nopre_GA
function [cdf_lf] = gradl_CDF_GA_ue(para_Pre, taugrid, nmixtures, Y, X)

    ncbeta = length(taugrid);
    [nocovar,nsample] = size(X);
    
    % extract c and beta 
    para = vec2mat(para_Pre(1:nocovar*ncbeta), ncbeta);
    
    % extract mu, lambda, sigma
    ParameterDist=para_Pre((nocovar*ncbeta+1):(nocovar*ncbeta+3*nmixtures-2));
    
    lambdapre=(ParameterDist([1:(nmixtures-1)]));
    mupre=(ParameterDist([(nmixtures):(2*nmixtures-2)]));
    sigma=ParameterDist([(2*nmixtures-1):end]);
    
    [lambda,mu,lambda3] = preprocesslambdamu(lambdapre,mupre);
    
    ParameterDist_nopre = [lambda, mu, sigma]; 
    
    [cdf_lf] = gradl_CDF_nopre_GA_ue(para, ParameterDist_nopre, taugrid, nmixtures, nocovar, Y, X);
    
    cdf_lf = -1*cdf_lf;
 
end
