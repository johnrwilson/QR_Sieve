function [lf,g] = gradllfnopreTdist(sv1,ntau, nsample, yt, x1r, x2r)
% Calculate the gradient of llf1
% To Calculate the gradient of llf1, it is necessary to calculate llf1. So
% might as well add llf1 as an output to this function.
%log-likelihood function, defined as llf1.

% Extract the parameters
beta0=sv1(1:ntau);
beta1=sv1((ntau+1):(2*ntau));
beta2=sv1((2*ntau+1):(3*ntau));
nu=sv1((3*ntau+1));

% Initialize the log likelihood
lf=0;
% Initialize the gradient for lambda, mu and sigma in the mixture of 3
% normal pdfs
gnu = zeros(1,1);

% Initialize the gradient for betas
gbeta0 = zeros(ntau,1);
gbeta1 = zeros(ntau,1);
gbeta2 = zeros(ntau,1);

erroreachsample = zeros(ntau,1);

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
        leachsample=leachsample+ 1/ntau*tpdf(erroreachsample(j_tau),nu);
    end
    lf = lf-log(leachsample)/nsample;
    
    
    % Now have the likelihood as variable leachsample. Need to calculate
    % the gradients.
    % 3) Calculate the gradient for nu
    gnueachsample = zeros(1,1);
    
    for j_tau = [1:ntau]
        u = erroreachsample(j_tau);
        gnueachsample = gnueachsample + 1/ntau*tpdf(u,nu)*(  (u^2*(nu+1))/(2*(u^2/nu+1)*nu^2)  -  log( (u^2/nu) + 1 )/2 );
        
        temp = (1 + u^2/nu)^(-(nu+1)/2);
        temp = temp/nu/pi/( (gamma(nu/2))^2 );
        
        temp1 = gammaderi( (nu+1)/2 )/2*sqrt( nu*pi )*gamma( nu/2 );
        temp1 = temp1 - gamma( (nu+1)/2 )*sqrt(pi)*( 1/(2*sqrt(nu))*gamma(nu/2)  + gammaderi( nu/2 )/2*sqrt(nu)   );
        
        gnueachsample = gnueachsample + 1/ntau*temp*temp1;
    end
    
    gnueachsample = gnueachsample/leachsample;
    gnu = gnu - gnueachsample/nsample;
    
    
    
    % 6) Calculate the gradient for beta1 and beta2
    gbeta0eachsample = zeros(ntau,1);
    gbeta1eachsample = zeros(ntau,1);
    gbeta2eachsample = zeros(ntau,1);
    for j_tau = [1:ntau]
        
        u = erroreachsample(j_tau);
        
        temp = gamma(  (nu+1)/2 );
        temp = temp/sqrt(nu*pi)/gamma( nu/2 );
        temp = temp*( -(nu+1)/2 );
        temp = temp*2*u/nu;
        temp = temp*( 1+u^2/nu )^( -(nu+3)/2 );
        
        gbeta0eachsample(j_tau) = temp;
        
        
        gbeta1eachsample(j_tau) = gbeta0eachsample(j_tau) ;
        gbeta2eachsample(j_tau) = gbeta0eachsample(j_tau) ;
        
        
    end
    gbeta0eachsample = gbeta0eachsample*(-1/ntau)*1/leachsample;
    gbeta1eachsample = gbeta1eachsample*(-1/ntau)*x1/leachsample;
    gbeta2eachsample = gbeta2eachsample*(-1/ntau)*x2/leachsample;
    
    gbeta0 = gbeta0 - gbeta0eachsample/nsample;
    gbeta1 = gbeta1 - gbeta1eachsample/nsample;
    gbeta2 = gbeta2 - gbeta2eachsample/nsample;
    
    
    
end
g = [gbeta0;gbeta1;gbeta2;gnu];

return;



