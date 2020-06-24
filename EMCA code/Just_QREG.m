% This codes doe
% Haoyang Liu
% 8/8/2018

clear
close all

% Part A) Definition of constants
taugrid_qreg = [0.1 0.25 0.5 0.75 0.9];
ntau = length(taugrid_qreg);
% Number of iterations (or runs of simulations)
iter = 500;
nsample = 100000;
%
n_epsilon = 3;
epsilon_v = [0,2,4];

recorder_qreg = nan(iter,n_epsilon,3*ntau);


% Part B) Estimation
for j_iter = [1:iter]
    j_iter
    
    % Betas are simulated in a continuous range. So no need to change that.
    tau_simu = rand(1,nsample);
    beta0_simu = tau_simu;
    beta1_simu = exp(tau_simu);
    beta2_simu = sqrt(tau_simu);
    
    x1r = exp(randn(1,nsample));
    x2r = exp(randn(1,nsample));
    
    y_n = normrnd(0,1,1,nsample);
    
    for j_sigma = [1:n_epsilon]
        y_s = beta0_simu + beta1_simu.*x1r + beta2_simu.*x2r;
        y = y_n*epsilon_v(j_sigma) + y_s;
        y = y';
        
        %% 3) QREG
        X = [ones(1,nsample); x1r; x2r]';
 
        
        [fit] = quantlsfVector(X,y,taugrid_qreg);
        fit_1 = [squeeze(fit(:,1));squeeze(fit(:,2));squeeze(fit(:,3))];
        recorder_qreg(j_iter,j_sigma,:)=fit_1';
        
    end
    
end

save just_qreg

bias = squeeze(mean(recorder_qreg,1)) - ones(3,1)*[taugrid_qreg, exp(taugrid_qreg), sqrt(taugrid_qreg)];

bias2 = bias(:,[ntau+1:2*ntau]);
bias3 = bias(:,[2*ntau+1:3*ntau]);

return;



