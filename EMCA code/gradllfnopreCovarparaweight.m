function [lf,g] = gradllfnopreCovarparaweight(sv1,ntau, nsample, nmixtures, yt, X, weight)
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
for j_sample = [1:nsample]
    
    xvector = X([j_sample],:);
    %%% x1 = x1r(j_sample);
    %%% x2 = x2r(j_sample);
    % 1) Calculate the error term in each sample corresponding to each ntau
    
    % 1*ntau
    erroreachsample = repmat(yt(j_sample),[1 ntau]) - xvector*betamatrix;
    
    
    % 2) Calculate the likelihood
    % Likelihood for each sample
    leachsample= (1/sum(weight))*sum( mnnoprevector(erroreachsample,ParameterDist,nmixtures).*weight );    
    lf = lf-log(leachsample)/nsample;
    
    
    % Now have the likelihood as variable leachsample. Need to calculate
    % the gradients.
    % 3) Calculate the gradient for lambda
    glambdaeachsample = zeros(nmixtures,1);
    for j_mixture = [1:nmixtures]
        glambdaeachsample(j_mixture) = glambdaeachsample(j_mixture) + (1/sum(weight))*sum(normpdf((erroreachsample-mu(j_mixture))/sigma(j_mixture)).*weight)/(sigma(j_mixture));
    end
    glambdaeachsample = glambdaeachsample/leachsample;
    glambda = glambda - glambdaeachsample/nsample;
    
    
    % 4) Calculate the gradient for mu
    gmueachsample = zeros(nmixtures,1);
    for j_mixture = [1:nmixtures]
        u_matrix(j_mixture,:) = (erroreachsample-mu(j_mixture))/sigma(j_mixture);
    end
    for j_mixture = [1:nmixtures]
        u_vector =  squeeze(u_matrix(j_mixture,:));
        gmueachsample(j_mixture) = gmueachsample(j_mixture) + (1/sum(weight))*sum(  lambda(j_mixture)*u_vector.*normpdf(u_vector)/((sigma(j_mixture))^2).*weight  );
    end
    gmueachsample = gmueachsample/leachsample;
    gmu = gmu - gmueachsample/nsample;
    
    
    % 5) Calculate the gradient for sigma
    gsigmaeachsample = zeros(nmixtures,1);
    for j_mixture = [1:nmixtures]
        
     u_vector =  squeeze(u_matrix(j_mixture,:));
        gsigmaeachsample(j_mixture) = gsigmaeachsample(j_mixture) - (1/sum(weight))*lambda(j_mixture)*sum( normpdf(u_vector).*weight )/((sigma(j_mixture))^2);
        gsigmaeachsample(j_mixture) = gsigmaeachsample(j_mixture) + (1/sum(weight))*lambda(j_mixture)*sum( (u_vector.^2).*normpdf(u_vector).*weight )/((sigma(j_mixture))^2);
        
    end
    gsigmaeachsample = gsigmaeachsample/leachsample;
    gsigma = gsigma - gsigmaeachsample/nsample;
    
    
    % 6) Calculate the gradient for beta1 and beta2
    gbeta0eachsample = zeros(ntau,1);
  
        for j_mixture = [1:nmixtures]

            u_vector =  squeeze(u_matrix(j_mixture,:));
   
            gbeta0eachsample = gbeta0eachsample  - (((lambda(j_mixture)*u_vector.*weight.*normpdf(u_vector)/((sigma(j_mixture))^2) ))/(sum(weight)))';
        end
        

    gbetamatrixeachsample =  repmat(gbeta0eachsample, [1 ncovar]);
    gbetamatrixeachsample = gbetamatrixeachsample.*repmat(xvector,[ntau 1])*(-1)/leachsample;
    
    
    gbetamatrix = gbetamatrix - gbetamatrixeachsample/nsample;
    
    
    
    
end
g = [reshape(gbetamatrix, ntau*ncovar, 1);glambda;gmu;gsigma];

return;



