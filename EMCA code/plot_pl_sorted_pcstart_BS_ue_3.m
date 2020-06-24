% This code plots the full sample piecewise constant and linear results with confidence interval from BS
% Real data
% Updated on 12/9/2018
clear;
close all;

%%% Part A) Set parameters and load data
year = '10';
switch year
    case '80'
        load RData_pl_Angrist80_sub1_ntau15_white1_ue1
        data=csvread('Angrist80.csv',1,0);
    case '90'
        load RData_pl_Angrist90_sub1_ntau15_white1_ue1
        data=csvread('Angrist90.csv',1,0);
    case '00'
        load RData_pl_Angrist00_sub1_ntau15_white1_ue1
         data=csvread('Angrist00.csv',1,0);
    otherwise
        load RData_pl_Jacob10_sub1_ntau15_white1_ue1
        load RData_pl_Jacob10_sub1_ntau15_white1_ue_ga
        data=csvread('Jacob10.csv',1,0);
end


[X_orig, y_orig] = white_subsample(data);
taugrid_qreg = [1:1:99]/100;
[fit] = quantlsfVector(X_orig,y_orig,taugrid_qreg);

beta_qreg_full = fit'; clear fit;
beta_WLS_start_sorted_full = beta_WLS_start_sorted; clear beta_WLS_start_sorted
beta_pc_s_full = beta_pc_s; clear beta_pc_s;

switch year
    case '80'
        load RData_pl_Angrist80_sub21_ntau15_white500_ue1
        postfixdata = 'Angrist 80'
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





% beta_qreg_std = squeeze(std(recorder_qreg_repeat_reshape));
beta_WLS_start_sorted_std = squeeze(std(recorder_WLS_start_repeat_sorted));
beta_pc_s_std = squeeze(std(recorder_pc_start_repeat_recon));

%beta_pc_s_std = squeeze(std(beta_pc_s_repeat));


% D-2) Plot
j_covar = 1;

figure;
hold on;

output_matrix_qreg = zeros(2, 99);
output_matrix_pc = zeros(4, 15);
output_matrix_pl = zeros(4, 15);
output_matrix_pl_1 = zeros(2,15);


plot(taugrid_qreg,beta_qreg_full(j_covar+1,:),'g','LineWidth',1);
output_matrix_qreg(1, :) = taugrid_qreg;
output_matrix_qreg(2, :) = beta_qreg_full(j_covar+1,:);

plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:),'k','LineWidth',1);
output_matrix_pc(1, :) = taugrid_midpoint;
output_matrix_pc(2, :) = beta_WLS_start_sorted_full(j_covar+1,:);

plot(taugrid_ue,beta_pc_s_full(j_covar+1,:),'r','LineWidth',1);
output_matrix_pl(1, :) = taugrid_ue;
output_matrix_pl(2, :) = beta_pc_s_full(j_covar+1,:);



output_matrix_pl_1(1,:) =  taugrid_midpoint;
pl_temp = beta_pc_s_full(j_covar+1,:);
pl_temp(1) = 2/3*(pl_temp(1)) + 1/3*(pl_temp(2));
pl_temp(end) = 2/3*(pl_temp(end)) + 1/3*(pl_temp(end-1));
output_matrix_pl_1(2,:) =  pl_temp;


plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:) + 1.96*beta_WLS_start_sorted_std(j_covar+1,:),'k--');
output_matrix_pc(3, :) = beta_WLS_start_sorted_full(j_covar+1,:) + 1.96*beta_WLS_start_sorted_std(j_covar+1,:);
plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:) - 1.96*beta_WLS_start_sorted_std(j_covar+1,:),'k--');
output_matrix_pc(4, :) = beta_WLS_start_sorted_full(j_covar+1,:) - 1.96*beta_WLS_start_sorted_std(j_covar+1,:);



% beta_pc_s_upper =  prctile(squeeze(recorder_pc_start_repeat_recon(:,j_covar+1,:)),97.5);
% beta_pc_s_lower =  prctile(squeeze(recorder_pc_start_repeat_recon(:,j_covar+1,:)),2.5);
    
    
plot(taugrid_ue,beta_pc_s_full(j_covar+1,:) + 1.96*beta_pc_s_std(j_covar+1,:),'r--');
output_matrix_pl(3, :) = beta_pc_s_full(j_covar+1,:) + 1.96*beta_pc_s_std(j_covar+1,:);
plot(taugrid_ue,beta_pc_s_full(j_covar+1,:) - 1.96*beta_pc_s_std(j_covar+1,:),'r--');
output_matrix_pl(4, :) = beta_pc_s_full(j_covar+1,:) - 1.96*beta_pc_s_std(j_covar+1,:);

% plot(taugrid_ue,beta_pc_s_upper,'m--');
% plot(taugrid_ue,beta_pc_s_lower,'m--');

legend('quantile regression','piecewise constant','piecewise linear(fmincon)');
title(sprintf('Data: %s. # of bootstrap samples: %d',postfixdata, n_repeatsample))
return;
output_matrix_qreg = [output_matrix_qreg]';
csvwrite(sprintf('%s qreg.csv',postfixdata),output_matrix_qreg,0,0);

output_matrix_pc = [output_matrix_pc]';
csvwrite(sprintf('%s pc.csv',postfixdata),output_matrix_pc,0,0);

output_matrix_pl = [output_matrix_pl]';
csvwrite(sprintf('%s pl.csv',postfixdata),output_matrix_pl,0,0);

output_matrix_pl_1 = [output_matrix_pl_1]';
csvwrite(sprintf('%s pl 1.csv',postfixdata),output_matrix_pl_1,0,0);


%axis([min(taugrid_midpoint) max(taugrid_midpoint) min(output_matrix(:, 6)) max(output_matrix(:, 5))])
print('-dpng','-r0',sprintf('BS_pl_ue_%s_%d_ext',postfix,n_repeatsample));









