function [cdf_lf, cdf_g] = gradl_CDF_nopre_1_ue(para, ParameterDist_nopre, taugrid, nmixtures, nocovar, N, Y, X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradl_CDF_nopre_1_ue
% Operationalize calculation of log likelihood function
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % A) Construct erroreach sample
    ntau = length(taugrid)-1;
    step_size_ue = taugrid([2:end]) - taugrid([1:end-1]);
    nsample = size(X,2);
    % Original para here is constants and slopes
    % Change it to constants and increments
    para_increments = para;
    para_increments(:,[2:end]) = para_increments(:,[2:end]).*(ones(nocovar,1)*step_size_ue);
    est_level = reconstruct_beta(para_increments);
    % Size of erroreachsample is # of observations times # of taus
    erroreachsample = repmat(Y, [(ntau+1) 1]) - est_level'*X;

    % B) Construct xB matrix (this section is clear and understood)
    % xB_k corresponds to xB_k the second to the last equation on page 1 of change_of_variables.pdf
    % Because it is the ratio in the change of variables procedure, it cannot be 0
    % xB_k = x0beta0^k+x1beta1^k+x2beta2^k+...
    % xB_k is slope multiplied by x, so its dimension is ntau * nsample 
    para_beta = para(:,2:end);
    xB_k = zeros(ntau, nsample);
    k=1;
    while k <= nocovar
        xB_k = xB_k + repmat(X(k,:),ntau,1).*repmat(para_beta(k,:)',1,nsample);
        k = k+1;
    end
 

    % C) Derivations independent of the taugrid
    f_u = 0*erroreachsample;
    F_u = 0*erroreachsample;
    F_u_r = 0*erroreachsample;
    G_u = zeros(ntau+1, nsample, nmixtures);
    G_u_r = zeros(ntau+1, nsample, nmixtures);
    H_u = zeros(ntau+1, nsample, nmixtures);
    P1_u = zeros(ntau+1, nsample, nmixtures);
    P1_u_r = zeros(ntau+1, nsample, nmixtures);
    P2_u = zeros(ntau+1, nsample, nmixtures);
    P2_u_r = zeros(ntau+1, nsample, nmixtures);
    C_u = 0*erroreachsample;
    
    lambda = ParameterDist_nopre(1:nmixtures);
    mu = ParameterDist_nopre(nmixtures+1: 2*nmixtures);
    sigma = ParameterDist_nopre(2*nmixtures+1: 3*nmixtures);
    
    % Variables that use the normcdf function:
    % F_u, G_u, and P1_u P2_u
    % All these variables were 

    for j = 1:nmixtures
        f_u = f_u+lambda(j)./sigma(j)*normpdf(erroreachsample,mu(j),sigma(j));
        F_u = F_u+lambda(j)*normcdf(erroreachsample,mu(j),sigma(j));
        % F_u_r
        F_u_r = F_u_r + lambda(j)*normcdf(-1*erroreachsample,mu(j),sigma(j)); 
        % G_u_r
        G_u(:,:,j) = normcdf(erroreachsample,mu(j),sigma(j));
        G_u_r(:,:,j) = normcdf(-1*erroreachsample,mu(j),sigma(j));
        H_u(:,:,j) = -lambda(j)./sigma(j).*normpdf((erroreachsample-mu(j))./sigma(j),0,1);
        
        P1_u(:,:,j)= -lambda(j)./sigma(j).*normcdf(erroreachsample,mu(j),sigma(j));
        P1_u_r(:,:,j)= -lambda(j)./sigma(j).*normcdf(-1*erroreachsample,mu(j),sigma(j));
        
        P2_u(:,:,j)= -lambda(j)./sigma(j).*(normpdf(erroreachsample,mu(j),sigma(j)).*(erroreachsample-mu(j))...
                    -normcdf(erroreachsample,mu(j),sigma(j))+1/2);
                
        P2_u_r(:,:,j)= -lambda(j)./sigma(j).*(normpdf(erroreachsample,mu(j),sigma(j)).*(erroreachsample-mu(j))...
                    +normcdf(-1*erroreachsample,mu(j),sigma(j))+1/2);     
                
        C_u = C_u-lambda(j)./sigma(j).*normpdf((erroreachsample-mu(j))./sigma(j),0,1);
    end
    
    f_u_diff = zeros(ntau, nsample);
    F_u_diff = zeros(ntau, nsample);
    F_u_diff_r = zeros(ntau, nsample);
    C_u_diff = zeros(ntau, nsample);
    
    for k = 1:ntau
        f_u_diff(k,:) = f_u(k,:) - f_u(k+1,:);
        F_u_diff(k,:) = F_u(k,:) - F_u(k+1,:);
        F_u_diff_r(k,:) = F_u_r(k,:) - F_u_r(k+1,:);
        C_u_diff(k,:) = C_u(k,:) - C_u(k+1,:);
    end
    
    
    F_u_diff = F_u_diff - ( F_u_diff == 0 ).* F_u_diff_r; 
   
    G_u_diff = zeros(ntau, nsample, nmixtures);
    G_u_diff_r = zeros(ntau, nsample, nmixtures);
    H_u_diff = zeros(ntau, nsample, nmixtures);
    P1_u_diff = zeros(ntau, nsample, nmixtures);
    P1_u_diff_r = zeros(ntau, nsample, nmixtures);
    P2_u_diff = zeros(ntau, nsample, nmixtures);
    P2_u_diff_r = zeros(ntau, nsample, nmixtures);
    for j = 1:nmixtures
        for k = 1:ntau
            G_u_diff(k,:,j) = G_u(k,:,j)-G_u(k+1,:,j);
            G_u_diff_r(k,:,j) = G_u_r(k,:,j)-G_u_r(k+1,:,j);
            H_u_diff(k,:,j) = H_u(k,:,j)-H_u(k+1,:,j);
            P1_u_diff(k,:,j) = P1_u(k,:,j)-P1_u(k+1,:,j);
            P1_u_diff_r(k,:,j) = P1_u_r(k,:,j)-P1_u_r(k+1,:,j);
            P2_u_diff(k,:,j) = P2_u(k,:,j)-P2_u(k+1,:,j);
            P2_u_diff_r(k,:,j) = P2_u_r(k,:,j)-P2_u_r(k+1,:,j);
        end
    end
    
    G_u_diff_temp = G_u_diff - ( G_u_diff == 0 ).* G_u_diff_r; 
    P1_u_diff = P1_u_diff - ( G_u_diff == 0 ).* P1_u_diff_r; 
    P2_u_diff = P2_u_diff + ( G_u_diff == 0 ).* P2_u_diff_r; 
    G_u_diff = G_u_diff_temp;
    
    
    F_value = F_u_diff./xB_k;
    G_value = G_u_diff./repmat(xB_k,1,1,3);
    H_value = H_u_diff./repmat(xB_k,1,1,3);
    P1_value = P1_u_diff./repmat(xB_k,1,1,3);
    P2_value = P2_u_diff./repmat(xB_k,1,1,3);
    C_value = C_u_diff./xB_k;
    
    cdf_lf_eachsample = sum(F_value);

    cdf_lf = sum(log(cdf_lf_eachsample))/nsample; %% log likelihood function
  
    cdf_glambdaeachsample = squeeze(sum(G_value, 1)./cdf_lf_eachsample);
    cdf_glambdaall = sum(cdf_glambdaeachsample)/nsample;

    cdf_gmueachsample = squeeze(sum(H_value, 1)./cdf_lf_eachsample);
    cdf_gmuall = sum(cdf_gmueachsample)/nsample;

    cdf_gsigmaeachsample = squeeze(sum(P1_value, 1)./cdf_lf_eachsample)+squeeze(sum(P2_value, 1)./cdf_lf_eachsample);
    cdf_gsigma = sum(cdf_gsigmaeachsample)/nsample;

    cdf_gceachsample = zeros(nmixtures, nsample);
    for k = 1:nocovar
        cdf_gceachsample(k,:) = sum(C_value)./cdf_lf_eachsample.*X(k,:);
    end
    cdf_gc = sum(cdf_gceachsample, 2)/nsample;



    % D) \partial log(likelihood) \partial beta. This one depends on the tau grid
    % Size of B1_value is ntau * nsample
    B1_value = (f_u(2:end,:).*xB_k.*(step_size_ue'*ones(1,nsample))-F_u_diff)./(xB_k.^2);
    B2_value = f_u_diff./xB_k;
    B2_value = cumsum(B2_value(2:end,:),'reverse').*(step_size_ue(1:end-1)'*ones(1,nsample));
    B2_value(ntau,:) = 0;
    cdf_gbetaeachsample = repmat((B1_value-B2_value)./cdf_lf_eachsample, 1,1,nocovar);
    for k = 1:nocovar
        cdf_gbetaeachsample(:,:,k)= cdf_gbetaeachsample(:,:,k).*X(k,:);
    end
    cdf_gbeta = (squeeze(sum(cdf_gbetaeachsample, 2))/nsample)';
    
    % E) Combine all partial derivatives
    % g_c_beta = [cdf_gc'; cdf_gbeta'];
    cdf_gbeta = cdf_gbeta';
    g_c_beta_v = [cdf_gc', cdf_gbeta(:)']; 
    cdf_g = [g_c_beta_v, cdf_glambdaall, cdf_gmuall, cdf_gsigma];
end
