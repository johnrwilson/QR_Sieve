clear 
close all

load Test_MLE_Residuals_Median_three_as3_logn_E
Para_true = [beta0_true, beta1_true, beta2_true, sv_true];

fval_recorder_truth = nan(iter,1);
fval_recorder_qreg_start_1 = nan(iter,1);
for j_iter = [1:iter]
    display('got here')

   
    x1r = x1matrix(j_iter,:);
    x2r = x2matrix(j_iter,:);
    y = ytmatrix(j_iter,:) ;        
    y = y';
  
    x = [ones(1,nsample);x1r;x2r]';
        
    fval_recorder_truth(j_iter) = gradllfCovarparavector(Para_true, ntau, nsample, nmixtures,1, y, x);
    fval_recorder_qreg_start_1(j_iter) = gradllfCovarparavector(squeeze(recorder_qreg_start(j_iter,:)), ntau, nsample, nmixtures,1, y, x);

end

y_all = reshape(ytmatrix,[iter*nsample,1]);
x1_all = reshape(x1matrix,[iter*nsample,1]);
x2_all = reshape(x2matrix,[iter*nsample,1]);

x_all = [ones(1,iter*nsample);x1_all';x2_all']';
fval_recorder_truth_all = gradllfCovarparavector(Para_true, ntau, iter*nsample, nmixtures,1, y_all, x_all);


fval_recorder_truth_qreg_start_all = nan(iter,1);
for j_iter = [1:iter]
    display('got here')
    fval_recorder_truth_qreg_start_all(j_iter) = gradllfCovarparavector(squeeze(recorder_qreg_start(j_iter,:)), ntau, iter*nsample, nmixtures,1, y_all, x_all);
end



save Test_MLE_Residuals_Median_three_as3_logn_E


return;


