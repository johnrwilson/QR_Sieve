function [lf,g] = gradllfnopre(sv1,ntau, nsample, nmixtures, yt, x1r, x2r)
    % Calculate the gradient of llf1
    % To Calculate the gradient of llf1, it is necessary to calculate llf1. So
    % might as well add llf1 as an output to this function.
    %log-likelihood function, defined as llf1.

    % Extract the parameters
    beta0=sv1(1:ntau);
    beta1=sv1((ntau+1):(2*ntau));
    beta2=sv1((2*ntau+1):(3*ntau));
    ParameterDist=sv1((3*ntau+1):(3*ntau+3*nmixtures));

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
    gbeta0 = zeros(ntau,1);
    gbeta1 = zeros(ntau,1);
    gbeta2 = zeros(ntau,1);

    % The main loop. Each iteration is for each sample. The formulas in
    % NoteMLE.pdf is for each iteration. Thus the gradient and llf1 are sum of
    % the formulas in NoteMLE over samples.
    % Don't forget to change the sign to negative
    for j_sample = [1:nsample]
        x1 = x1r(j_sample);
        x2 = x2r(j_sample);
        % 1) Calculate the error term in each sample corresponding to each ntau
        for j_tau = [1:ntau]
            erroreachsample(j_tau) = (yt(j_sample)-x1.*beta1(j_tau)-x2.*beta2(j_tau)-beta0(j_tau));
        end

        % 2) Calculate the likelihood
        % Likelihood for each sample
        leachsample=0;
        % \int_{0}^1 f_{\epsilon}( y - x \beta(\ntau)) d \ntau,
        for j_tau = [1:ntau]
            leachsample=leachsample+ 1/ntau*mnnopre(erroreachsample(j_tau),ParameterDist,nmixtures);
        end
        lf = lf-log(leachsample)/nsample;


        % Now have the likelihood as variable leachsample. Need to calculate
        % the gradients.
        % 3) Calculate the gradient for lambda
        glambdaeachsample = zeros(nmixtures,1);
        for j_mixture = [1:nmixtures]
            for j_tau = [1:ntau]
                glambdaeachsample(j_mixture) = glambdaeachsample(j_mixture) + 1/ntau*normpdf((erroreachsample(j_tau)-mu(j_mixture))/sigma(j_mixture))/(sigma(j_mixture));
            end
        end
        glambdaeachsample = glambdaeachsample/leachsample;
        glambda = glambda - glambdaeachsample/nsample;

        % 4) Calculate the gradient for mu
        gmueachsample = zeros(nmixtures,1);
        for j_mixture = [1:nmixtures]
            for j_tau = [1:ntau]
                u = (erroreachsample(j_tau)-mu(j_mixture))/sigma(j_mixture);
                gmueachsample(j_mixture) = gmueachsample(j_mixture) + 1/ntau*lambda(j_mixture)*u*normpdf(u)/((sigma(j_mixture))^2);
            end
        end
        gmueachsample = gmueachsample/leachsample;
        gmu = gmu - gmueachsample/nsample;


        % 5) Calculate the gradient for sigma
        gsigmaeachsample = zeros(nmixtures,1);
        for j_mixture = [1:nmixtures]
            for j_tau = [1:ntau]
                u = (erroreachsample(j_tau)-mu(j_mixture))/sigma(j_mixture);
                gsigmaeachsample(j_mixture) = gsigmaeachsample(j_mixture) - 1/ntau*lambda(j_mixture)*normpdf(u)/((sigma(j_mixture))^2);
                gsigmaeachsample(j_mixture) = gsigmaeachsample(j_mixture) + 1/ntau*lambda(j_mixture)*u^2*normpdf(u)/((sigma(j_mixture))^2);
            end
        end
        gsigmaeachsample = gsigmaeachsample/leachsample;
        gsigma = gsigma - gsigmaeachsample/nsample;


        % 6) Calculate the gradient for beta1 and beta2
       gbeta0eachsample = zeros(ntau,1);
       gbeta1eachsample = zeros(ntau,1);
       gbeta2eachsample = zeros(ntau,1);
        for j_tau = [1:ntau]
            for j_mixture = [1:nmixtures]
                u = (erroreachsample(j_tau)-mu(j_mixture))/sigma(j_mixture);
                gbeta0eachsample(j_tau) = gbeta0eachsample(j_tau)  - lambda(j_mixture)*u*normpdf(u)/((sigma(j_mixture))^2);
                gbeta1eachsample(j_tau) = gbeta1eachsample(j_tau)  - lambda(j_mixture)*u*normpdf(u)/((sigma(j_mixture))^2);
                gbeta2eachsample(j_tau) = gbeta2eachsample(j_tau)  - lambda(j_mixture)*u*normpdf(u)/((sigma(j_mixture))^2);
            end

        end
        gbeta0eachsample = gbeta0eachsample*(-1/ntau)*1/leachsample;        
        gbeta1eachsample = gbeta1eachsample*(-1/ntau)*x1/leachsample;    
        gbeta2eachsample = gbeta2eachsample*(-1/ntau)*x2/leachsample;    

        gbeta0 = gbeta0 - gbeta0eachsample/nsample;  
        gbeta1 = gbeta1 - gbeta1eachsample/nsample;
        gbeta2 = gbeta2 - gbeta2eachsample/nsample;    



    end
    g = [gbeta0;gbeta1;gbeta2;glambda;gmu;gsigma];

return;



