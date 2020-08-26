function [lf,g] = gradllfCovarparavector(parametervectorPre, ntau, nmixtures, nslice, yt, X)
% Calculate the gradient of llf1
% To Calculate the gradient of llf1, it is necessary to calculate llf1. So
% might as well add llf1 as an output to this function.

[nsample,ncovar] = size(X);
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

[ begin_index,end_index,index_length ] = dividematrix( nsample,nslice );
if nslice >1
    parfor j_slice = [1:nslice]
        index_v = ([begin_index(j_slice):end_index(j_slice) ]);
        [lf_v(j_slice),g_v(:,j_slice)] = gradllfnopreCovarparavector(parametervector,ntau, index_length(j_slice), nmixtures, yt(index_v), X(index_v,:));
    end
else
    j_slice = 1;
    index_v = ([begin_index(j_slice):end_index(j_slice) ]);
    [lf_v(j_slice),g_v(:,j_slice)] = gradllfnopreCovarparavector(parametervector,ntau, index_length(j_slice), nmixtures, yt(index_v), X(index_v,:));
    
end
lf = sum(lf_v.*index_length)/nsample;
g = g_v*(index_length')/nsample;

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

g = [g([1:ntau*ncovar]);glambda;gmu;gsigma]';

return;



