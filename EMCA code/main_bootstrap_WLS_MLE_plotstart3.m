clear;
close all;

year = 1980;
postfixload = 'Angrist80_33';

[taugrid,BS_start,WLS_start,qreg_start,qreg,BS_start_smooth,WLS_start_smooth] = loadsmoothbootstrap(year);

load MLE_WLS_Vector50_Angrist80_33_sub19 recorder_WLS_start_repeat 
load MLE_real_data_qregAngrist80mixturethree_33 X

BS = mean(recorder_WLS_start_repeat,1);

X_mean = mean(X);
BS_reshape = reshape(BS(1:ntau*ncovar),ntau,ncovar);
[BS_sort] = sortbeta_1(X_mean,BS_sort,ntau,ncovar);

return;

hold on;
plot(taugrid,qreg,'k');
plot(taugrid,qreg,'k');

%plot(taugrid,betahatsmoothBS(:,5),'k--');
