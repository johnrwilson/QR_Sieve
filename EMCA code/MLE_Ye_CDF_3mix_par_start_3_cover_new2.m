% This code is the second code calculating the coverage probability
% It constructs the lower and the upper bounds on those parameters
% Haoyang Liu
% 11/15/2018

clear;
close all;

load Ye_GA_3mix_3_par_start_3_cover_new_5

taugrid_true = [0:(ntau-1)]/ntau + 1/(2*(ntau-1));
%taugrid_true = [0:(ntau-1)]/(ntau-1);
beta_true = [b0(taugrid_true),b1(taugrid_true),b2(taugrid_true)];


for j_iter = [1:iter]
    recorder_qreg_start_sorted_temp = squeeze(recorder_qreg_start_bootstrap_sorted(j_iter,:,:));
    recorder_lambda_sorted_temp =  squeeze(recorder_lambda_sorted_bootstrap(j_iter,:,:));
    recorder_mu_sorted_temp = squeeze(recorder_mu_sorted_bootstrap(j_iter,:,:));
    recorder_sigma_sorted_temp = squeeze(recorder_sigma_sorted_bootstrap(j_iter,:,:));
    
    std_beta(j_iter,:) = std(recorder_qreg_start_sorted_temp,1);
    recorder_qreg_start_sorted_lower(j_iter,:) = -1.96*std(recorder_qreg_start_sorted_temp,1) + recorder_qreg_start_sorted(j_iter, :);
    recorder_qreg_start_sorted_upper(j_iter,:) = 1.96*std(recorder_qreg_start_sorted_temp,1) + recorder_qreg_start_sorted(j_iter, :);
    
    recorder_lambda_sorted_lower(j_iter,:) = -1.96*std(recorder_lambda_sorted_temp,1) + recorder_lambda_sorted(j_iter, :);
    recorder_lambda_sorted_upper(j_iter,:) = 1.96*std(recorder_lambda_sorted_temp,1) + recorder_lambda_sorted(j_iter, :);

        
    recorder_mu_sorted_lower(j_iter,:) = -1.96*std(recorder_mu_sorted_temp,1) + recorder_mu_sorted(j_iter, :);
    recorder_mu_sorted_upper(j_iter,:) = 1.96*std(recorder_mu_sorted_temp,1) + recorder_mu_sorted(j_iter, :);
    
    recorder_sigma_sorted_lower(j_iter,:) = -1.96*std(recorder_sigma_sorted_temp,1) + recorder_sigma_sorted(j_iter, :);
    recorder_sigma_sorted_upper(j_iter,:) = 1.96*std(recorder_sigma_sorted_temp,1) + recorder_sigma_sorted(j_iter, :);
 
    
end

cp_beta = mean((recorder_qreg_start_sorted_lower<=beta_true).*(recorder_qreg_start_sorted_upper>=beta_true));
mean(cp_beta)

cp_lambda = mean((recorder_lambda_sorted_lower<=lambdavector_true).*(recorder_lambda_sorted_upper>=lambdavector_true));
mean(cp_lambda)

cp_mu = mean((recorder_mu_sorted_lower<=muvector_true).*(recorder_mu_sorted_upper>=muvector_true));
mean(cp_mu)

cp_sigma =  mean((recorder_sigma_sorted_lower<=[1,1,1]).*(recorder_sigma_sorted_upper>=[1,1,1]));
mean(cp_sigma)


mean([cp_beta, cp_lambda, cp_mu, cp_sigma])


