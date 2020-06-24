function [lf,g] = gradllfnopreCovarparavector_1(sv1,ntau, nsample, nmixtures, yt, X)
% Calculate the gradient of llf1
% To Calculate the gradient of llf1, it is necessary to calculate llf1. So
% might as well add llf1 as an output to this function.
%log-likelihood function, defined as llf1.
[temp,ncovar] = size(X);

betamatrix = vec2mat(sv1([1:(ntau*ncovar)]),ntau);
% Extract the parameters
%%%  beta0=sv1(1:ntau);
%%%  beta1=sv1((ntau+1):(2*ntau));
%%%  beta2=sv1((2*ntau+1):(3*ntau));


ParameterDist=sv1((ntau*ncovar+1):(ntau*ncovar+3*nmixtures));

% Preprocessed lambda. Lambda=weights of components, only need to specify
% the first nmixtures-1 weights.
lambda=(ParameterDist([1:(nmixtures)]));
% Preprocessed mu. mu= mean of components, only need to specify the first nmixtures-1 means.
mu=(ParameterDist([(nmixtures+1):(2*nmixtures)]));
%sigma= st.d of each component, need to specify for all nmixtures components.
sigma=ParameterDist([(2*nmixtures+1):(3*nmixtures)]);

% Initialize the log likelihood
lf=0;
% Initialize the gradient for lambda, mu and sigma in the mixture of 3
% normal pdfs
glambda = zeros(nmixtures,1);
gmu = zeros(nmixtures,1);
gsigma = zeros(nmixtures,1);
% Initialize the gradient for betas
gbetamatrix = zeros(ntau,ncovar);
gbeta0 = zeros(ntau,1);
%%%   gbeta1 = zeros(ntau,1);
%%%   gbeta2 = zeros(ntau,1);

% The main loop. Each iteration is for each sample. The formulas in
% NoteMLE.pdf is for each iteration. Thus the gradient and llf1 are sum of
% the formulas in NoteMLE over samples.
% Don't forget to change the sign to negative





    

    
    erroreachsample = repmat(yt,[1 ntau]) - X*betamatrix;
    
    
    % 2) Calculate the likelihood
    % Likelihood for each sample
    leachsample= 1/ntau*sum( mnnoprevectorvector(erroreachsample,ParameterDist,nmixtures),2);    
    lf =  - sum(log(leachsample))/nsample;
    
    
    % Now have the likelihood as variable leachsample. Need to calculate
    % the gradients.
    % 3) Calculate the gradient for lambda
    glambdaeachsample = zeros(nmixtures,nsample);
    for j_mixture = [1:nmixtures]
        glambdaeachsample(j_mixture,:) = glambdaeachsample(j_mixture,:) + (1/ntau*sum(normpdf((erroreachsample-mu(j_mixture))/sigma(j_mixture)),2)/(sigma(j_mixture)))';
    end
    glambdaeachsample = glambdaeachsample./(repmat(leachsample',[nmixtures 1]));
    glambda = - sum(glambdaeachsample,2)/nsample;
    
    
    % 4) Calculate the gradient for mu
    gmueachsample = zeros(nmixtures,nsample);
    for j_mixture = [1:nmixtures]
        u_matrix(j_mixture,:,:) = (erroreachsample-mu(j_mixture))/sigma(j_mixture);
    end
    for j_mixture = [1:nmixtures]
        u_vector =  squeeze(u_matrix(j_mixture,:,:));
        gmueachsample(j_mixture,:) = gmueachsample(j_mixture,:) + 1/ntau*(sum(  lambda(j_mixture)*u_vector.*normpdf(u_vector),2)/((sigma(j_mixture))^2)  )';
    end
    gmueachsample = gmueachsample./repmat(leachsample',[nmixtures 1]);
    gmu = - sum(gmueachsample,2)/nsample;
    
    
    % 5) Calculate the gradient for sigma
    gsigmaeachsample = zeros(nmixtures,nsample);
    for j_mixture = [1:nmixtures]
        
     u_vector =  squeeze(u_matrix(j_mixture,:,:));
        gsigmaeachsample(j_mixture,:) = gsigmaeachsample(j_mixture,:) - 1/ntau*lambda(j_mixture)*(sum( normpdf(u_vector),2 )/((sigma(j_mixture))^2))';
        gsigmaeachsample(j_mixture,:) = gsigmaeachsample(j_mixture,:) + 1/ntau*lambda(j_mixture)*(sum( (u_vector.^2).*normpdf(u_vector),2 )/((sigma(j_mixture))^2))';
        
    end
    gsigmaeachsample = gsigmaeachsample./repmat(leachsample',[nmixtures 1]);
    gsigma = - sum(gsigmaeachsample,2)/nsample;
    
    
    % 6) Calculate the gradient for beta1 and beta2
    gbeta0eachsample = zeros(nsample,ntau);
  
        for j_mixture = [1:nmixtures]
            u_vector =  squeeze(u_matrix(j_mixture,:,:));
            gbeta0eachsample = gbeta0eachsample  - ( lambda(j_mixture)*u_vector.*normpdf(u_vector)/((sigma(j_mixture))^2) );
        end
        

    gbetamatrixeachsample = zeros(nsample,ntau,ncovar);
        
    for j_covar = [1:ncovar]
        gbetamatrixeachsample(:,:,j_covar) = gbeta0eachsample;
    end
        

    gbetamatrixeachsample = permute(gbetamatrixeachsample,[1,3,2]);
    
    xtemp = zeros(nsample,ncovar,ntau);
    ltemp = xtemp;
    for j_tau = [1:ntau]
        xtemp(:,:,j_tau) = X;
        for j_covar = [1:ncovar]
            ltemp(:,j_covar,j_tau) = leachsample;
        end
    end
        
    gbetamatrixeachsample = gbetamatrixeachsample.*xtemp*(-1/ntau)./ltemp;
    
    
    gbetamatrix = - squeeze(sum(gbetamatrixeachsample,1))/nsample;
    
    
   
g = [reshape(gbetamatrix', ntau*ncovar, 1);glambda;gmu;gsigma];

return;



