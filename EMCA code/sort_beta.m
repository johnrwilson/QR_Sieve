function [sorted_beta] = sort_beta(X_mean,fit_hat,ncovar,ntau) 

    beta_qreg_start =  reshape(fit_hat(1,[1 : ncovar*ntau])',ntau,ncovar);
    fval_tau_qreg_start = X_mean*beta_qreg_start';
    [V,I_qreg_start] = sort(fval_tau_qreg_start);
    for j_x = [1:3]
        sorted_beta(1,[(j_x-1)*ntau+1 : j_x*ntau]) = beta_qreg_start(I_qreg_start,j_x)';
    end
    
end
    
    
    
    
