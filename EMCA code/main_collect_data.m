clear;
close all;

%WLS

load MLE_real_data_WLSAngrist80mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 5;
X_mean = mean(X);

recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);


[num_obs,num_x] = size(X);
num_obs_WLS_80 = num_obs;
fval_recorder_WLS_start_80 = fval_recorder_WLS_start;
beta_WLS_start_sorted_80 = beta_WLS_start_sorted;
exit_recorder_WLS_start_80 = exit_recorder_WLS_start;


load MLE_real_data_WLSAngrist90mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 5;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_WLS_90 = num_obs;
fval_recorder_WLS_start_90 = fval_recorder_WLS_start;
beta_WLS_start_sorted_90 = beta_WLS_start_sorted;
exit_recorder_WLS_start_90 = exit_recorder_WLS_start;

load MLE_real_data_WLSAngrist00mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 5;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_WLS_00 = num_obs;
fval_recorder_WLS_start_00 = fval_recorder_WLS_start;
beta_WLS_start_sorted_00 = beta_WLS_start_sorted;
exit_recorder_WLS_start_00 = exit_recorder_WLS_start;

load MLE_real_data_WLSJacob10mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 5;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_WLS_10 = num_obs;
fval_recorder_WLS_start_10 = fval_recorder_WLS_start;
beta_WLS_start_sorted_10 = beta_WLS_start_sorted;
exit_recorder_WLS_start_10 = exit_recorder_WLS_start;

%qreg

load MLE_real_data_qregAngrist80mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 5;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_80 = num_obs;
fval_recorder_qreg_start_80 = fval_recorder_qreg_start;
beta_qreg_start_sorted_80 = beta_qreg_start_sorted;
exit_recorder_qreg_start_80 = exit_recorder_qreg_start;


load MLE_real_data_qregAngrist90mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 5;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_90 = num_obs;
fval_recorder_qreg_start_90 = fval_recorder_qreg_start;
beta_qreg_start_sorted_90 = beta_qreg_start_sorted;
exit_recorder_qreg_start_90 = exit_recorder_qreg_start;



load MLE_real_data_qregAngrist00mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 5;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);


[num_obs,num_x] = size(X);
num_obs_qreg_00 = num_obs;
fval_recorder_qreg_start_00 = fval_recorder_qreg_start;
beta_qreg_start_sorted_00 = beta_qreg_start_sorted;
exit_recorder_qreg_start_00 = exit_recorder_qreg_start;


load MLE_real_data_qregJacob10mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 5;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);


[num_obs,num_x] = size(X);
num_obs_qreg_10 = num_obs;
fval_recorder_qreg_start_10 = fval_recorder_qreg_start;
beta_qreg_start_sorted_10 = beta_qreg_start_sorted;
exit_recorder_qreg_start_10 = exit_recorder_qreg_start;


%WLS White


load MLE_real_data_WLSAngristWhite80mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 4;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);


[num_obs,num_x] = size(X);
num_obs_WLS_80White = num_obs;
fval_recorder_WLS_start_80White = fval_recorder_WLS_start;
beta_WLS_start_sorted_80White = beta_WLS_start_sorted;
exit_recorder_WLS_start_80White = exit_recorder_WLS_start;


load MLE_real_data_WLSAngristWhite90mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 4;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);


[num_obs,num_x] = size(X);
num_obs_WLS_90White = num_obs;
fval_recorder_WLS_start_90White = fval_recorder_WLS_start;
beta_WLS_start_sorted_90White = beta_WLS_start_sorted;
exit_recorder_WLS_start_90White = exit_recorder_WLS_start;

load MLE_real_data_WLSAngristWhite00mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 4;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_WLS_00White = num_obs;
fval_recorder_WLS_start_00White = fval_recorder_WLS_start;
beta_WLS_start_sorted_00White = beta_WLS_start_sorted;
exit_recorder_WLS_start_00White = exit_recorder_WLS_start;


load MLE_real_data_WLSJacobWhite10mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 4;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_WLS_10White = num_obs;
fval_recorder_WLS_start_10White = fval_recorder_WLS_start;
beta_WLS_start_sorted_10White = beta_WLS_start_sorted;
exit_recorder_WLS_start_10White = exit_recorder_WLS_start;

% qreg White

load MLE_real_data_qregAngristWhite80mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 4;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_80White = num_obs;
fval_recorder_qreg_start_80White = fval_recorder_qreg_start;
beta_qreg_start_sorted_80White = beta_qreg_start_sorted;
exit_recorder_qreg_start_80White = exit_recorder_qreg_start;




load MLE_real_data_qregAngristWhite90mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 4;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_90White = num_obs;
fval_recorder_qreg_start_90White = fval_recorder_qreg_start;
beta_qreg_start_sorted_90White = beta_qreg_start_sorted;
exit_recorder_qreg_start_90White = exit_recorder_qreg_start;


load MLE_real_data_qregAngristWhite00mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 4;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_00White = num_obs;
fval_recorder_qreg_start_00White = fval_recorder_qreg_start;
beta_qreg_start_sorted_00White = beta_qreg_start_sorted;
exit_recorder_qreg_start_00White = exit_recorder_qreg_start;


load MLE_real_data_qregJacobWhite10mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 4;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_10White = num_obs;
fval_recorder_qreg_start_10White = fval_recorder_qreg_start;
beta_qreg_start_sorted_10White = beta_qreg_start_sorted;
exit_recorder_qreg_start_10White = exit_recorder_qreg_start;



%WLS Black


load MLE_real_data_WLSAngristBlack80mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 4;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_WLS_80Black = num_obs;
fval_recorder_WLS_start_80Black = fval_recorder_WLS_start;
beta_WLS_start_sorted_80Black = beta_WLS_start_sorted;
exit_recorder_WLS_start_80Black = exit_recorder_WLS_start;


load MLE_real_data_WLSAngristBlack90mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 4;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_WLS_90Black = num_obs;
fval_recorder_WLS_start_90Black = fval_recorder_WLS_start;
beta_WLS_start_sorted_90Black = beta_WLS_start_sorted;
exit_recorder_WLS_start_90Black = exit_recorder_WLS_start;

load MLE_real_data_WLSAngristBlack00mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 4;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);


[num_obs,num_x] = size(X);
num_obs_WLS_00Black = num_obs;
fval_recorder_WLS_start_00Black = fval_recorder_WLS_start;
beta_WLS_start_sorted_00Black = beta_WLS_start_sorted;
exit_recorder_WLS_start_00Black = exit_recorder_WLS_start;


load MLE_real_data_WLSJacobBlack10mixturethree_33 recorder_WLS_start exit_recorder_WLS_start ntau X fval_recorder_WLS_start;
nx = 4;
X_mean = mean(X);
recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*nx),ntau,nx);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,nx);


[num_obs,num_x] = size(X);
num_obs_WLS_10Black = num_obs;
fval_recorder_WLS_start_10Black = fval_recorder_WLS_start;
beta_WLS_start_sorted_10Black = beta_WLS_start_sorted;
exit_recorder_WLS_start_10Black = exit_recorder_WLS_start;

% qreg Black


load MLE_real_data_qregAngristBlack80mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 4;
X_mean = mean(X);

recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_80Black = num_obs;
fval_recorder_qreg_start_80Black = fval_recorder_qreg_start;
beta_qreg_start_sorted_80Black = beta_qreg_start_sorted;
exit_recorder_qreg_start_80Black = exit_recorder_qreg_start;




load MLE_real_data_qregAngristBlack90mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 4;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_90Black = num_obs;
fval_recorder_qreg_start_90Black = fval_recorder_qreg_start;
beta_qreg_start_sorted_90Black = beta_qreg_start_sorted;
exit_recorder_qreg_start_90Black = exit_recorder_qreg_start;


load MLE_real_data_qregAngristBlack00mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 4;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);

[num_obs,num_x] = size(X);
num_obs_qreg_00Black = num_obs;
fval_recorder_qreg_start_00Black = fval_recorder_qreg_start;
beta_qreg_start_sorted_00Black = beta_qreg_start_sorted;
exit_recorder_qreg_start_00Black = exit_recorder_qreg_start;


load MLE_real_data_qregJacobBlack10mixturethree_33 recorder_qreg_start exit_recorder_qreg_start ntau X fval_recorder_qreg_start;
nx = 4;
X_mean = mean(X);
recorder_qreg_start_reshape = reshape(recorder_qreg_start(1:ntau*nx),ntau,nx);
[beta_qreg_start_sorted] = sortbeta_1(X_mean,recorder_qreg_start_reshape,ntau,nx);


[num_obs,num_x] = size(X);
num_obs_qreg_10Black = num_obs;
fval_recorder_qreg_start_10Black = fval_recorder_qreg_start;
beta_qreg_start_sorted_10Black = beta_qreg_start_sorted;
exit_recorder_qreg_start_10Black = exit_recorder_qreg_start;


fval_recorder_qreg(1) = fval_recorder_qreg_start_80;
fval_recorder_qreg(2) = fval_recorder_qreg_start_90;
fval_recorder_qreg(3) = fval_recorder_qreg_start_00;
fval_recorder_qreg(4) = fval_recorder_qreg_start_10;

fval_recorder_qreg(5) = fval_recorder_qreg_start_80White;
fval_recorder_qreg(6) = fval_recorder_qreg_start_90White;
fval_recorder_qreg(7) = fval_recorder_qreg_start_00White;
fval_recorder_qreg(8) = fval_recorder_qreg_start_10White;

fval_recorder_qreg(9) = fval_recorder_qreg_start_80Black;
fval_recorder_qreg(10) = fval_recorder_qreg_start_90Black;
fval_recorder_qreg(11) = fval_recorder_qreg_start_00Black;
fval_recorder_qreg(12) = fval_recorder_qreg_start_10Black;

fval_recorder_WLS(1) = fval_recorder_WLS_start_80;
fval_recorder_WLS(2) = fval_recorder_WLS_start_90;
fval_recorder_WLS(3) = fval_recorder_WLS_start_00;
fval_recorder_WLS(4) = fval_recorder_WLS_start_10;

fval_recorder_WLS(5) = fval_recorder_WLS_start_80White;
fval_recorder_WLS(6) = fval_recorder_WLS_start_90White;
fval_recorder_WLS(7) = fval_recorder_WLS_start_00White;
fval_recorder_WLS(8) = fval_recorder_WLS_start_10White;

fval_recorder_WLS(9) = fval_recorder_WLS_start_80Black;
fval_recorder_WLS(10) = fval_recorder_WLS_start_90Black;
fval_recorder_WLS(11) = fval_recorder_WLS_start_00Black;
fval_recorder_WLS(12) = fval_recorder_WLS_start_10Black;



num_obs_qreg(1) = num_obs_qreg_80;
num_obs_qreg(2) = num_obs_qreg_90;
num_obs_qreg(3) = num_obs_qreg_00;
num_obs_qreg(4) = num_obs_qreg_10;

num_obs_qreg(5) = num_obs_qreg_80White;
num_obs_qreg(6) = num_obs_qreg_90White;
num_obs_qreg(7) = num_obs_qreg_00White;
num_obs_qreg(8) = num_obs_qreg_10White;


num_obs_qreg(9) = num_obs_qreg_80Black;
num_obs_qreg(10) = num_obs_qreg_90Black;
num_obs_qreg(11) = num_obs_qreg_00Black;
num_obs_qreg(12) = num_obs_qreg_10Black;

num_obs_WLS(1) = num_obs_WLS_80;
num_obs_WLS(2) = num_obs_WLS_90;
num_obs_WLS(3) = num_obs_WLS_00;
num_obs_WLS(4) = num_obs_WLS_10;

num_obs_WLS(5) = num_obs_WLS_80White;
num_obs_WLS(6) = num_obs_WLS_90White;
num_obs_WLS(7) = num_obs_WLS_00White;
num_obs_WLS(8) = num_obs_WLS_10White;


num_obs_WLS(9) = num_obs_WLS_80Black;
num_obs_WLS(10) = num_obs_WLS_90Black;
num_obs_WLS(11) = num_obs_WLS_00Black;
num_obs_WLS(12) = num_obs_WLS_10Black;






coeffcients(:,1) = squeeze(beta_WLS_start_sorted_80(:,2))'
coeffcients(:,2) = squeeze(beta_WLS_start_sorted_90(:,2))'
coeffcients(:,3) = squeeze(beta_WLS_start_sorted_00(:,2))'
coeffcients(:,4) = squeeze(beta_WLS_start_sorted_10(:,2))'
coeffcients(:,5) = squeeze(beta_qreg_start_sorted_80(:,2))'
coeffcients(:,6) = squeeze(beta_qreg_start_sorted_90(:,2))'
coeffcients(:,7) = squeeze(beta_qreg_start_sorted_00(:,2))'
coeffcients(:,8) = squeeze(beta_qreg_start_sorted_10(:,2))'



coeffcients(:,9) = squeeze(beta_WLS_start_sorted_80White(:,2))'
coeffcients(:,10) = squeeze(beta_WLS_start_sorted_90White(:,2))'
coeffcients(:,11) = squeeze(beta_WLS_start_sorted_00White(:,2))'
coeffcients(:,12) = squeeze(beta_WLS_start_sorted_10White(:,2))'
coeffcients(:,13) = squeeze(beta_qreg_start_sorted_80White(:,2))'
coeffcients(:,14) = squeeze(beta_qreg_start_sorted_90White(:,2))'
coeffcients(:,15) = squeeze(beta_qreg_start_sorted_00White(:,2))'
coeffcients(:,16) = squeeze(beta_qreg_start_sorted_10White(:,2))'


coeffcients(:,17) = squeeze(beta_WLS_start_sorted_80Black(:,2))'
coeffcients(:,18) = squeeze(beta_WLS_start_sorted_90Black(:,2))'
coeffcients(:,19) = squeeze(beta_WLS_start_sorted_00Black(:,2))'
coeffcients(:,20) = squeeze(beta_WLS_start_sorted_10Black(:,2))'
coeffcients(:,21) = squeeze(beta_qreg_start_sorted_80Black(:,2))'
coeffcients(:,22) = squeeze(beta_qreg_start_sorted_90Black(:,2))'
coeffcients(:,23) = squeeze(beta_qreg_start_sorted_00Black(:,2))'
coeffcients(:,24) = squeeze(beta_qreg_start_sorted_10Black(:,2))'








