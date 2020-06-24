% This code implements the new gradient function based on the codees by Lei
% para_Pre is the parameter vector prior to inferrring the omitted error
% mixture
% Other input variables are named in a straight forward way. 
function [lf,g] = gradl_Lei_orig(para_Pre, taugrid, nmixtures, yt, X)

% Note that ntau is the number of segments plus 1. 
% The number of parameters for each covar is also the number of segments plus 1 
ncbeta = length(taugrid)+1;
[nocovar,nsample] = size(X);

% Extract the parameters
ParameterDist=para_Pre((nocovar*ncbeta+1):(nocovar*ncbeta+3*nmixtures-2));

% Preprocessed lambda. Lambda=weights of components, only need to specify
% the first nmixtures-1 weights.
lambdapre=(ParameterDist([1:(nmixtures-1)]));
% Preprocessed mu. mu= mean of components, only need to specify the first nmixtures-1 means.
mupre=(ParameterDist([(nmixtures):(2*nmixtures-2)]));
%sigma= st.d of each component, need to specify for all nmixtures components.
sigma=ParameterDist([(2*nmixtures-1):end]);

% Prerocess lambda and mu
[lambda,mu,lambda3] = preprocesslambdamu(lambdapre,mupre);



para = vec2mat(para_Pre(1:ncbeta*nocovar), (ncbeta));
 %%%%%%% Don't forget to change this!!!!, The input of gradl_nopre_Lei has
 %%%%%%% been changed.
 

ParameterDist_nopre = [lambda, mu, sigma]; 
[lf, g] = gradl_nopre_Lei_orig(para, ParameterDist_nopre, taugrid, nmixtures, nocovar, nsample, yt, X);
 
glambdaall = g([(nocovar*ncbeta+1):(nocovar*ncbeta+nmixtures)])';
gmuall     = g([(nocovar*ncbeta+nmixtures+1):(nocovar*ncbeta+2*nmixtures)])';
gsigma  = g([(nocovar*ncbeta+2*nmixtures+1):(nocovar*ncbeta+3*nmixtures)])';
if nmixtures>1
    glambda = glambdaall([1:(nmixtures-1)]) - glambdaall(nmixtures)-( mupre'*lambda3 + (mupre*(lambdapre')))/(lambda3^2)*gmuall(nmixtures) ;
    gmu =  gmuall([1:(nmixtures-1)]) - lambdapre'/(lambda3)*gmuall(nmixtures);
end
if nmixtures==1
    glambda = [];
    gmu =[];
end

lf = -1*lf;
g = -1*[g([1:ncbeta*nocovar])'; glambda; gmu; gsigma]';

return;






