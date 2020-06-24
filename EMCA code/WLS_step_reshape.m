function [recorder_WLS,recorder_WLS_sort] = WLS_step(fit_qreg,X,y,n_WLS_iter)

[ntau,ncovar] = size(fit_qreg);
nsample = length(y);
Y = repmat(y, [1 ntau]);

% Preallocation for weight matrix
weightedX = nan(nsample,ncovar);
betahat = nan(ncovar,ntau);


ehat = Y-X*(fit_qreg');
sigmahat=std(y-X*quantlsf(X,y,.5)); % use only std of median regression residuals



for j_WLS=1:n_WLS_iter
    pdfmatrix = normpdf(ehat./sigmahat);
    % the weight matrix is also nsample * ntau
    weight=ntau*pdfmatrix./repmat(sum(pdfmatrix,2),[1 ntau]);

    for j_tau = [1:ntau]
        weightV = weight(:,j_tau)';
        for j_covar = [1:ncovar]
            weightedX(:,j_covar) = weightV.*(X(:,j_covar)');
        end
        weightedy = weightV'.*y;
        betahat(:,j_tau)=(X'*weightedX)\(X'*weightedy);
    end
    
    ehat=Y-X*betahat;
    sigmahat=(sum(sum(weight.*(ehat.^2)))/nsample/ntau)^.5;
    
end
recorder_WLS = reshape(betahat',ncovar*ntau,1)';
recorder_WLS_sort = sortbeta_1(mean(X),recorder_WLS,ntau,ncovar);


