% This code plots the full sample piecewise constant and linear results with confidence interval from BS
% Real data
% Updated on 12/9/2018
clear;
close all;

%%% Part A) Set parameters and load data
year = '90';
switch year
    case '80'
        load RData_pl_Angrist80_sub1_ntau15_white1_ue1
    case '90'
        load RData_pl_Angrist90_sub1_ntau15_white1_ue1
    case '00'
        load RData_pl_Angrist00_sub1_ntau15_white1_ue1
    otherwise
        load RData_pl_Jacob10_sub1_ntau15_white1_ue1
end


beta_qreg_full = fit'; clear fit;
beta_WLS_start_sorted_full = beta_WLS_start_sorted; clear beta_WLS_start_sorted
beta_pc_s_full = beta_pc_s; clear beta_pc_s;

switch year
    case '80'
        load RData_pl_Angrist80_sub21_ntau15_white500_ue1
        postfixdata =  'Angrist 80'
    case '90'
        load RData_pl_Angrist90_sub21_ntau15_white500_ue1
        postfixdata =  'Angrist 90'
    case '00'
        load RData_pl_Angrist00_sub21_ntau15_white500_ue1
        postfixdata =  'Angrist 00'
    otherwise
        load RData_pl_Jacob10_sub21_ntau15_white500_ue1
        postfixdata =  'Jacob 10'
        
end





beta_qreg_std = squeeze(std(recorder_qreg_repeat_reshape));
beta_WLS_start_sorted_std = squeeze(std(recorder_WLS_start_repeat_sorted));
beta_pc_s_std = squeeze(std(recorder_pc_start_repeat_recon));

%beta_pc_s_std = squeeze(std(beta_pc_s_repeat));


% D-2) Plot
j_covar = 1;

figure;
hold on;

output_matrix = zeros(6, 15);

plot(taugrid,beta_qreg_full(j_covar+1,:),'g','LineWidth',1);
output_matrix(1, :) = beta_qreg_full(j_covar+1,:);
output_matrix(2, :) = beta_qreg_full(j_covar+1,:);
% plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:),'k','LineWidth',1);
% output_matrix(, :) = beta_WLS_start_sorted_full(j_covar+1,:);
plot(taugrid_ue,beta_pc_s_full(j_covar+1,:),'r','LineWidth',1);
output_matrix(3, :) = taugrid_ue;
output_matrix(4, :) = beta_pc_s_full(j_covar+1,:);


% plot(taugrid,beta_qreg_full(j_covar+1,:) + 1.96*beta_qreg_std(j_covar+1,:),'g--');
% output_matrix(, :) = beta_qreg_full(j_covar+1,:) + 1.96*beta_qreg_std(j_covar+1,:);
% plot(taugrid,beta_qreg_full(j_covar+1,:) - 1.96*beta_qreg_std(j_covar+1,:),'g--');
%output_matrix(, :) = beta_qreg_full(j_covar+1,:) - 1.96*beta_qreg_std(j_covar+1,:);

% plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:) + 1.96*beta_WLS_start_sorted_std(j_covar+1,:),'k--');
% output_matrix(, :) = beta_WLS_start_sorted_full(j_covar+1,:) + 1.96*beta_WLS_start_sorted_std(j_covar+1,:);
% plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:) - 1.96*beta_WLS_start_sorted_std(j_covar+1,:),'k--');
% output_matrix(, :) = beta_WLS_start_sorted_full(j_covar+1,:) - 1.96*beta_WLS_start_sorted_std(j_covar+1,:);



beta_pc_s_upper =  prctile(squeeze(recorder_pc_start_repeat_recon(:,j_covar+1,:)),97.5);
beta_pc_s_lower =  prctile(squeeze(recorder_pc_start_repeat_recon(:,j_covar+1,:)),2.5);
    
    
plot(taugrid_ue,beta_pc_s_full(j_covar+1,:) + 1.96*beta_pc_s_std(j_covar+1,:),'r--');
output_matrix(5, :) = beta_pc_s_full(j_covar+1,:) + 1.96*beta_pc_s_std(j_covar+1,:);
plot(taugrid_ue,beta_pc_s_full(j_covar+1,:) - 1.96*beta_pc_s_std(j_covar+1,:),'r--');
output_matrix(6, :) = beta_pc_s_full(j_covar+1,:) - 1.96*beta_pc_s_std(j_covar+1,:);

plot(taugrid_ue,beta_pc_s_upper,'m--');
plot(taugrid_ue,beta_pc_s_lower,'m--');

legend('quantile regression','piecewise constant','piecewise linear(fmincon)');
title(sprintf('Data: %s. # of bootstrap samples: %d',postfixdata, n_repeatsample))

output_matrix = [output_matrix]';
csvwrite(sprintf('%s_1.csv',postfixdata),output_matrix,1,0);

axis([min(taugrid) max(taugrid) min(output_matrix(:, 6)) max(output_matrix(:, 5))])

print('-dpng','-r0',sprintf('BS_pl_ue_%s_%d_1',postfix,n_repeatsample));









