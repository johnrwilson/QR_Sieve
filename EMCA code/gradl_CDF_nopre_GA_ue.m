% Function gradl_CDF_nopre_GA_ue
% Almost the same as gradl_CDF_nopre_GA, but can deal with 
function [cdf_lf] = gradl_CDF_nopre_GA_ue(para, ParameterDist_nopre, taugrid, nmixtures, nocovar, Y, X)

    ntau = length(taugrid) - 1;
    step_size_ue = taugrid([2:end]) - taugrid([1:end-1]);
    
    nsample = size(X,2);

    % This reconstruct_beta function converts increments to the values at
    % those knots
    beta = reconstruct_beta(para);

    % Size of erroreachsample is # of observations times # of taus
    erroreachsample = repmat(Y, [(ntau+1) 1]) - beta'*X;
    
    % xB_k = x0beta0^k+x1beta1^k+x2beta2^k+...
    para_beta = para(:,2:end)./(ones(nocovar,1)*step_size_ue);
    xB_k = zeros(ntau, nsample);
   
    for k = 1:nocovar
        xB_k = xB_k + repmat(X(k,:),ntau,1).*repmat(para_beta(k,:)',1,nsample);
    end
    
    F_u = 0*erroreachsample;
    F_u_r = 0*erroreachsample;
   
    lambda = ParameterDist_nopre(1:nmixtures);
    mu = ParameterDist_nopre(nmixtures+1: 2*nmixtures);
    sigma = ParameterDist_nopre(2*nmixtures+1: 3*nmixtures);
    
    % Variables that use the normcdf function:

    for j = 1:nmixtures      
        F_u = F_u + lambda(j)*normcdf(erroreachsample,mu(j),sigma(j));
        % F_u_r
        F_u_r = F_u_r + lambda(j)*normcdf(-1*erroreachsample,mu(j),sigma(j)); 
    end  
   
    F_u_diff = zeros(ntau, nsample);
    F_u_diff_r = zeros(ntau, nsample);
    
    for k = 1:ntau    
        F_u_diff(k,:) = F_u(k,:) - F_u(k+1,:);
        F_u_diff_r(k,:) = F_u_r(k,:) - F_u_r(k+1,:);
    end

    
    F_value = F_u_diff./xB_k;
   
    cdf_lf_eachsample = sum(F_value);

    cdf_lf = sum(log(cdf_lf_eachsample + 0.0000001))/nsample; %% log likelihood function
  
end
