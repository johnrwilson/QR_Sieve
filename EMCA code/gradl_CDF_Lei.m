%% 4) Function gradl_CDF_Lei
function [cdf_lf,cdf_g] = gradl_CDF_Lei(para_Pre, taugrid, nmixtures, Y, X)

    ncbeta = length(taugrid)+1;
    [nocovar,nsample] = size(X);
    
    % extract c and beta 
    para = vec2mat(para_Pre(nocovar+1:nocovar*ncbeta), ncbeta-1);
    c = para_Pre(1:nocovar)';
    para = [c para];
    
    % extract mu, lambda, sigma
    ParameterDist=para_Pre((nocovar*ncbeta+1):(nocovar*ncbeta+3*nmixtures-2));
    
    lambdapre=(ParameterDist([1:(nmixtures-1)]));
    mupre=(ParameterDist([(nmixtures):(2*nmixtures-2)]));
    sigma=ParameterDist([(2*nmixtures-1):end]);
    
    [lambda,mu,lambda3] = preprocesslambdamu(lambdapre,mupre);
    
    ParameterDist_nopre = [lambda, mu, sigma]; 
    
    [cdf_lf, cdf_g] = gradl_CDF_nopre_Lei_1(para, ParameterDist_nopre, taugrid, nmixtures, nocovar, nsample, Y, X);
    
    glambdaall = cdf_g([(nocovar*ncbeta+1):(nocovar*ncbeta+nmixtures)])';
    gmuall     = cdf_g([(nocovar*ncbeta+nmixtures+1):(nocovar*ncbeta+2*nmixtures)])';
    gsigma  = cdf_g([(nocovar*ncbeta+2*nmixtures+1):(nocovar*ncbeta+3*nmixtures)])';
    if nmixtures>1
        glambda = glambdaall([1:(nmixtures-1)]) - glambdaall(nmixtures)-( mupre'*lambda3 + (mupre*(lambdapre')))/(lambda3^2)*gmuall(nmixtures) ;
        gmu =  gmuall([1:(nmixtures-1)]) - lambdapre'/(lambda3)*gmuall(nmixtures);
    end
    if nmixtures==1
        glambda = [];
        gmu =[];
    end

    cdf_lf = -1*cdf_lf;
    cdf_g = -1*[cdf_g([1:ncbeta*nocovar])'; glambda; gmu; gsigma]';
end
