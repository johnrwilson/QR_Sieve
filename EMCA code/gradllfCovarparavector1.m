function [lf,g] = gradllfCovarparavector1(parametervectorPre,ntau, nsample, nmixtures, nslice, yt,X)
% Calculate the gradient of llf1
% To Calculate the gradient of llf1, it is necessary to calculate llf1. So
% might as well add llf1 as an output to this function.

[temp,ncovar] = size(X);
%log-likelihood function, defined as llf1.

% Extract the parameters

ParameterDist=parametervectorPre((ncovar*ntau+1):(ncovar*ntau+3*nmixtures-1));

% Preprocessed lambda. Lambda=weights of components, only need to specify
% the first nmixtures-1 weights.
lambda=(ParameterDist([1:(nmixtures)]));
% Preprocessed mu. mu= mean of components, only need to specify the first nmixtures-1 means.
mupre=(ParameterDist([(nmixtures+1):(2*nmixtures-1)]));
%sigma= st.d of each component, need to specify for all nmixtures components.
sigma=ParameterDist([(2*nmixtures):end]);

% Prerocess lambda and mu
[mu] = preprocesslambdamu1(lambda,mupre);

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
    
    glambdafirst = glambdaall([1:(nmixtures-1)]) - mupre'./lambda(nmixtures).*gmuall(nmixtures);
    glambdalast = glambdaall(nmixtures) - mu(nmixtures)/lambda(nmixtures)*gmuall(nmixtures);
    glambda = [glambdafirst;glambdalast];
    
    
      gmu =  gmuall([1:(nmixtures-1)]) - lambda([1:(nmixtures-1)])'/(lambda(nmixtures))*gmuall(nmixtures);
end



if nmixtures==1
    glambda = glambdaall;
    gmu =[];
end

g = [g([1:ntau*ncovar]);glambda;gmu;gsigma];

return;



