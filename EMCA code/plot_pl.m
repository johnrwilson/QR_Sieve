% This code plots results from 


clear

close all 

load MLE_pw_linear_truth

taugrid_qreg = [0:1/9:1];
% Part A) True coefficients
beta1_true = b1(taugrid_qreg);
beta2_true = b2(taugrid_qreg);

% Part B) Fill in the parametes for qreg and start of qreg
beta1_qreg = recorder_qreg(:,[11:20]);
beta2_qreg = recorder_qreg(:,[21:30]);

beta_c1_qreg_start = recorder_qreg_start(:,[2, 13:21]);
beta_c2_qreg_start = recorder_qreg_start(:,[3, 22:30]);


beta1_qreg_start = pl_pc(beta_c1_qreg_start, 1/9);
beta2_qreg_start = pl_pc(beta_c2_qreg_start, 1/9);


bias_beta1_qreg = mean(beta1_qreg) - beta1_true;
bias_beta1_qreg_start = mean(beta1_qreg_start) - beta1_true;


bias_beta2_qreg = mean(beta2_qreg) - beta2_true;
bias_beta2_qreg_start = mean(beta2_qreg_start) - beta2_true;

figure; hold on;

plot(taugrid_qreg,bias_beta1_qreg,'b');

plot(taugrid_qreg,bias_beta1_qreg_start,'r');


legend('Quantile Regression','Piecewise Linear MLE after truth'); 
title('Average Bias in Beta 1 Estimates')
print('-dpng','-r0','beta1_pl_truth');


figure; hold on;

plot(taugrid_qreg,bias_beta2_qreg,'b');

plot(taugrid_qreg,bias_beta2_qreg_start,'r');


legend('Quantile Regression','Piecewise Linear MLE after truth'); 
title('Average Bias in Beta 2 Estimates')
print('-dpng','-r0','beta2_pl_truth');

