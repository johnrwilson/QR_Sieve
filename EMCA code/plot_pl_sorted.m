% This is the plotting code
% It plots sorted estimates
% It doesn't let sorting done at a finer grid 
% 1/23/2018

clear
close all 

%% 1) Choose among different results
truth_start = 0;
switch truth_start
    case -4 
        load pl_CDF_truth_c2_33_new
        print_string = "Truth except C_1";
    case -3
        load pl_CDF_qreg_3_33
        print_string = "QREG with True W2";
    case -2
        load pl_CDF_qreg_para
        print_string = "Only True Dist";
    case -1
        load pl_CDF_truth_para
        print_string = "Only True Beta";
    case 1
        load pl_CDF_truth
        print_string = "Truth";
    case 2
        load pl_unc_truth_c2_33
        print_string = "Truth except c_1, fminunc";
    case 3
        load pl_unc_truth_33
        print_string = "Truth, fminunc";
    case 4
        load pl_unc_qreg_para_33
        print_string = "qreg and true EIV parameters, fminunc";   
    
    
    otherwise
        %load pl_CDF_qreg_33_new
        load pl_CDF_qreg_33_1norm
        print_string = "QREG 33";
end

%% 2) Create "taugrid_break_sorted", which repeats each middle entry twice 
step_original = taugrid_qreg(3) - taugrid_qreg(2);
% taugrid_break is still a little different from taugrid_qreg because the
% first entry is 0 instead of 0+0.00001
taugrid_break = [0:step_original:1];
taugrid_break_sorted = repeat_HL(taugrid_break,2);
% Repeating it twice gives us end points
taugrid_break_sorted = taugrid_break_sorted(1:end-1);

%% 3) Keep replicates with valid results
recorder_qreg = recorder_qreg(~isnan(recorder_qreg_start(:,1)),:);
recorder_qreg_start = recorder_qreg_start(~isnan(recorder_qreg_start(:,1)),:);
[iter, temp] = size(recorder_qreg_start);

%% 4) Sort the estimates
% Part A) True coefficients
beta1_true = b1(taugrid_qreg);
beta2_true = b2(taugrid_qreg);

% Part B) Fill in the parametes for qreg and start of qreg
beta1_qreg = recorder_qreg(:,[(ntau+2):(2*ntau+2)]);
beta2_qreg = recorder_qreg(:,[(2*ntau+3):(3*ntau+3)]);

beta_c0_qreg_start = recorder_qreg_start(:,[1, 4:(ntau+3)]);
beta_c1_qreg_start = recorder_qreg_start(:,[2, (ntau+4):(2*ntau+3)]);
beta_c2_qreg_start = recorder_qreg_start(:,[3, (2*ntau+4):(3*ntau+3)]);

beta0_qreg_start = pl_pc(beta_c0_qreg_start, step_original);
beta1_qreg_start = pl_pc(beta_c1_qreg_start, step_original);
beta2_qreg_start = pl_pc(beta_c2_qreg_start, step_original);

% The following four variables are the average estimates before sorting
beta1_qreg_mean = mean(beta1_qreg);
beta1_qreg_start_mean = mean(beta1_qreg_start);
beta2_qreg_mean = mean(beta2_qreg);
beta2_qreg_start_mean = mean(beta2_qreg_start);

% C) Sorting 
beta0_qreg_start_mid = (beta0_qreg_start(:, [1:end-1]) + beta0_qreg_start(:, [1:end-1]))/2;
beta1_qreg_start_mid = (beta1_qreg_start(:, [1:end-1]) + beta1_qreg_start(:, [1:end-1]))/2;
beta2_qreg_start_mid = (beta2_qreg_start(:, [1:end-1]) + beta2_qreg_start(:, [1:end-1]))/2;

beta0_qreg_start_sorted = nan(iter, floor(2/step_original+0.00001));
beta1_qreg_start_sorted = nan(iter, floor(2/step_original+0.00001));
beta2_qreg_start_sorted = nan(iter, floor(2/step_original+0.00001));

for j_iter = [1:iter]
    x1r = squeeze(x1matrix(j_iter,:));
    x2r = squeeze(x2matrix(j_iter,:));
    
    X_mean = [1,mean(x1r),mean(x2r)];
     
    fval_tau_qreg_start = X_mean*[squeeze(beta0_qreg_start_mid(j_iter,:)); squeeze(beta1_qreg_start_mid(j_iter,:)); squeeze(beta2_qreg_start_mid(j_iter,:))];
    [V,I_qreg_start(j_iter,:)] = sort(fval_tau_qreg_start);
    
    beta0_qreg_start_sorted(j_iter,[1:2:end]) = beta0_qreg_start(j_iter,I_qreg_start(j_iter,:));
    beta0_qreg_start_sorted(j_iter,[2:2:end]) = beta0_qreg_start(j_iter,I_qreg_start(j_iter,:)+1);
    
    beta1_qreg_start_sorted(j_iter,[1:2:end]) = beta1_qreg_start(j_iter,I_qreg_start(j_iter,:));
    beta1_qreg_start_sorted(j_iter,[2:2:end]) = beta1_qreg_start(j_iter,I_qreg_start(j_iter,:)+1);
        
    beta2_qreg_start_sorted(j_iter,[1:2:end]) = beta2_qreg_start(j_iter,I_qreg_start(j_iter,:));
    beta2_qreg_start_sorted(j_iter,[2:2:end]) = beta2_qreg_start(j_iter,I_qreg_start(j_iter,:)+1);
    
end

beta1_qreg_start_sorted_mean = mean(beta1_qreg_start_sorted);
beta2_qreg_start_sorted_mean = mean(beta2_qreg_start_sorted);

%% 4) Plot
figure; hold on;


plot(taugrid_qreg,beta1_true,'g');
plot(taugrid_qreg,beta1_qreg_mean,'b');

plot(taugrid_break,beta1_qreg_start_mean,'r');
plot(taugrid_break_sorted,beta1_qreg_start_sorted_mean,'r--');

legend('Truth', 'Quantile Regression','Piecewise Linear MLE', 'Piecewise Linear MLE (Sorted)'); 
title(sprintf('Average Beta 1 Estimates (Start: %s)',print_string))
print('-dpng','-r0',sprintf('beta1_pl_%s',print_string));


figure; hold on;
plot(taugrid_qreg,beta2_true,'g');
plot(taugrid_qreg,beta2_qreg_mean,'b');

plot(taugrid_break,beta2_qreg_start_mean,'r');
plot(taugrid_break_sorted,beta2_qreg_start_sorted_mean,'r--');

legend('Truth', 'Quantile Regression','Piecewise Linear MLE', 'Piecewise Linear MLE (Sorted)'); 
title(sprintf('Average Beta 2 Estimates (Start: %s)',print_string))

print('-dpng','-r0',sprintf('beta2_pl_%s',print_string));



