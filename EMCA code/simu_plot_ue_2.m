
clear
close all


% Dummy for if a case involves piecewise linear
pl = 0;

% EIV = '3mix 99 knot PC';
% EIV = "T PC"
EIV = "Laplace PC"
% A) Load data
switch EIV
    case '3mix w WLS'
        load Ye_GA_3mix_3_WLS_ue_combined
        load Ye_GA_3mix_3_WLS_ue_99knot_500 taugrid recorder_qreg_repeat_reshape
        pl = 1
    case '2mix PC'
        % Estimating a mixture 2 normals as a mixture of 3 normals.
        load Ye_GA_2mix_3_WLS_ue
        load Ye_GA_2mix_3_WLS_ue_qreg taugrid recorder_qreg_repeat_reshape
    case 'Laplace PC'
        % Estimating a mixture 2 normals as a mixture of 3 normals.
        load Ye_GA_L_3_WLS_ue_600
        recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted([1:500], :, :);
        load Ye_GA_L_3_WLS_ue_qreg_500 taugrid recorder_qreg_repeat_reshape
    case '3mix dummy PC'
        load Ye_GA_3mix_3_WLS_ue_dummy
        load Ye_GA_3mix_3_WLS_ue_dummy_qreg taugrid recorder_qreg_repeat_reshape
    case 'T PC'
        load Ye_GA_T_3_WLS_ue_combined
        recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted(exit_recorder_WLS_start_repeat ==  1, :, :);
        recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted([1:500], :, :);
        load Ye_GA_T_3_WLS_ue_qreg_500 taugrid recorder_qreg_repeat_reshape
    case '3mix 99 knot PC'
        load Ye_GA_3mix_3_WLS_ue_99knot_500
        pl = 2;
  
    otherwise
end

taugrid_test= [.1:.1:.9];
beta_true = [b0(taugrid_test);b1(taugrid_test);b2(taugrid_test)];

if pl == 1
    recorder_MLE = recorder_pc_start_repeat_recon;
else
    recorder_MLE = recorder_WLS_start_repeat_sorted;
end
recorder_qreg_test = recorder_qreg_repeat_reshape(:,:,[10:10:90]);

recorder_MLE_test = nan(min(500,n_repeatsample),3,9);

for j = [1:2:9]
    recorder_MLE_test(:,:,j) = recorder_MLE(:,:,2+(j-1)/2*3);
end

for j = [2:2:8]
    recorder_MLE_test(:,:,j) = (1/2)*(recorder_MLE(:,:,3+(j-2)/2*3)+recorder_MLE(:,:,4+(j-2)/2*3));
end

if pl == 2
    recorder_MLE_test = recorder_MLE(:,:,[10:10:90]);
end

mean_qreg = squeeze(mean((recorder_qreg_test - permute(repmat(beta_true,[1 1 500]),[3 1 2]))));
mean_MLE = squeeze(mean((recorder_MLE_test - permute(repmat(beta_true,[1 1 min(n_repeatsample, 500) ]),[3 1 2]))));
output_matrix = [mean_qreg;mean_MLE]';
output_matrix = output_matrix(:,[1,4,2,5,3,6]);
csvwrite(sprintf('mean_%s.csv',EIV), output_matrix);

MSE_qreg = squeeze(mean((recorder_qreg_test - permute(repmat(beta_true,[1 1 500 ]),[3 1 2])).^2));
MSE_MLE = squeeze(mean((recorder_MLE_test - permute(repmat(beta_true,[1 1 min(n_repeatsample, 500) ]),[3 1 2])).^2));
output_matrix1 = [MSE_qreg;MSE_MLE]';
output_matrix1 = output_matrix1(:,[1,4,2,5,3,6]);
csvwrite(sprintf('mse_%s.csv',EIV), output_matrix1);
