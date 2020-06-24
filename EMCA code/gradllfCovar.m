function [lf,g] = gradllfCovar(parametervectorPre,ntau, nsample, nmixtures, yt,X)
% Calculate the gradient of llf1
% To Calculate the gradient of llf1, it is necessary to calculate llf1. So
% might as well add llf1 as an output to this function.

    [temp,ncovar] = size(X);
%log-likelihood function, defined as llf1.

% Extract the parameters

ParameterDist=parametervectorPre((ncovar*ntau+1):(ncovar*ntau+3*nmixtures-2));

% Preprocessed lambda. Lambda=weights of components, only need to specify
% the first nmixtures-1 weights.
lambdapre=(ParameterDist([1:(nmixtures-1)]));
% Preprocessed mu. mu= mean of components, only need to specify the first nmixtures-1 means.
mupre=(ParameterDist([(nmixtures):(2*nmixtures-2)]));
%sigma= st.d of each component, need to specify for all nmixtures components.
sigma=ParameterDist([(2*nmixtures-1):end]);

% Prerocess lambda and mu
[lambda,mu,lambda3] = preprocesslambdamu(lambdapre,mupre);

parametervector = [parametervectorPre(1:ntau*ncovar),lambda,mu,sigma];

[lf,g] = gradllfnopreCovar(parametervector,ntau, nsample, nmixtures, yt, X);


glambdaall = g([(ncovar*ntau+1):(ncovar*ntau+nmixtures)]);
gmuall     = g([(ncovar*ntau+nmixtures+1):(ncovar*ntau+2*nmixtures)]);
gsigma  = g([(ncovar*ntau+2*nmixtures+1):(ncovar*ntau+3*nmixtures)]);
if nmixtures>1
    glambda = glambdaall([1:(nmixtures-1)]) - glambdaall(nmixtures)-( mupre'*lambda3 + (mupre*(lambdapre')))/(lambda3^2)*gmuall(nmixtures) ;  
    gmu =  gmuall([1:(nmixtures-1)]) - lambdapre'/(lambda3)*gmuall(nmixtures);
end
if nmixtures==1
    glambda = [];
    gmu =[];
end

g = [g([1:ntau*ncovar]);glambda;gmu;gsigma];

return;



