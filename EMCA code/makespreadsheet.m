clear all;
close all;
ntau = 99;

% The iteration to calculate standard errors
j_iter = 30;
j_dataset = 1;
datamatrix = nan(ntau*3*3,4);
V_all = nan(5,5,ntau,9);
for j_pool = [1:3]
    
    if j_pool == 3
        postfix = 'Angrist';
    end
    if j_pool == 1
        postfix = 'AngristWhite';
    end
    if j_pool == 2
        postfix = 'AngristBlack';
    end
    
    for j_year = [1:3]
        if j_year == 1
            year = '80';
        end
        if j_year == 2
            year = '90';
        end
        if j_year == 3
            year = '00';
        end
        
        % Load the corresponding data set.
        load(sprintf('WLS_60_%s%s_99',postfix,year),'betahatMatrix','X','y','ncovar','sigmahat_v');
        
        % Get the beta coefficients
        X_mean = mean(X);
        [beta_WLS_sorted,I_WLS_sorted] = sortbeta_1(X_mean,squeeze(betahatMatrix(:,:,j_iter))',ntau,ncovar);
        
        % Calculate the row indeces in the data matrix
        j_begin = 1 + (j_dataset-1)*ntau;
        j_end = j_dataset*ntau;
        
        % Put the betas at the right place
        datamatrix([j_begin:j_end],1) = squeeze(beta_WLS_sorted(:,2));
        
        
        Y = repmat(y, [1 ntau]);
        
        betahat_previous =  betahatMatrix(:,:,(j_iter-1));
        sigmahat_previous = sigmahat_v(j_iter-1);
        ehat_previous = Y-X*betahat_previous;
        % the weight matrix is also nsample * ntau
        weight=ntau*normpdf(ehat_previous./sigmahat_previous)./repmat(sum(normpdf(ehat_previous./sigmahat_previous),2),[1 ntau]);
        
        betahat = betahatMatrix(:,:,j_iter);
        ehat=Y-X*betahat;
        sigmahat = sigmahat_v(j_iter);
        
        V = nan(ncovar,ncovar,ntau);
        for j_tau = [1:ntau]
            Wk = sparse(diag(weight(:,j_tau)));
            Ehat_V = sparse(diag(ehat(:,j_tau).^2));
            V_temp =  (X'*Wk*X)^(-1);
            V(:,:,j_tau) = V_temp*X'*Wk*Ehat_V*Wk*X*V_temp;
        end
        V_all([1:ncovar],[1:ncovar],:,j_dataset) = V;
        
        datamatrix([j_begin:j_end],2) = squeeze(V(2,2,I_WLS_sorted).^(1/2))';
        j_dataset = j_dataset+1;
    end
end

clear Wk Ehat_V betahatMatrix X y ncovar sigmahat_v ehat ehat_previous weight
save summarystd

return;


