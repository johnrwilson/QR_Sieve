clear


% The case for piecewise constant for a T EIV
% Load the first dataset 
load Ye_GA_T_3_WLS_ue
recorder_WLS_start_repeat_sorted_temp = recorder_WLS_start_repeat_sorted;
n_repeatsample_temp = n_repeatsample;
exit_recorder_WLS_start_repeat_temp = exit_recorder_WLS_start_repeat;

% Load the second dataset 
load Ye_GA_T_3_WLS_ue_1 recorder_WLS_start_repeat_sorted n_repeatsample exit_recorder_WLS_start_repeat;
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
exit_recorder_WLS_start_repeat_temp(1,[j_begin:j_start]) = exit_recorder_WLS_start_repeat;

% Load the nth dataset 
load Ye_GA_T_3_WLS_ue_2 recorder_WLS_start_repeat_sorted n_repeatsample exit_recorder_WLS_start_repeat;
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
exit_recorder_WLS_start_repeat_temp(1,[j_begin:j_start]) = exit_recorder_WLS_start_repeat;

% Load the nth dataset 
load Ye_GA_T_3_WLS_ue_3 recorder_WLS_start_repeat_sorted n_repeatsample exit_recorder_WLS_start_repeat;
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
exit_recorder_WLS_start_repeat_temp(1,[j_begin:j_start]) = exit_recorder_WLS_start_repeat;

% Load the nth dataset 
load Ye_GA_T_3_WLS_ue_4 recorder_WLS_start_repeat_sorted n_repeatsample exit_recorder_WLS_start_repeat;
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
exit_recorder_WLS_start_repeat_temp(1,[j_begin:j_start]) = exit_recorder_WLS_start_repeat;

% Load the nth dataset 
load Ye_GA_T_3_WLS_ue_5 recorder_WLS_start_repeat_sorted n_repeatsample exit_recorder_WLS_start_repeat;
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
exit_recorder_WLS_start_repeat_temp(1,[j_begin:j_start]) = exit_recorder_WLS_start_repeat;

clear n_repeatsample exit_recorder_WLS_start_repeat recorder_WLS_start_repeat_sorted

n_repeatsample = n_repeatsample_temp;
exit_recorder_WLS_start_repeat = exit_recorder_WLS_start_repeat_temp;
recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted_temp;

save Ye_GA_T_3_WLS_ue_combined






% The case for piecewise linear estimation for a Laplace EIV
% Load tau grid and ncovar
load Ye_GA_3mix_3_WLS_ue_3 taugrid_ue taugrid taugrid_midpoint ncovar

% Load the first set of data
load Ye_GA_3mix_3_WLS_ue_3 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
recorder_qreg_repeat_reshape_temp = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_3 recorder_pc_start_repeat
recorder_pc_start_repeat_temp = recorder_pc_start_repeat;
n_repeatsample_temp = n_repeatsample;

% Load the second set of data
load Ye_GA_3mix_3_WLS_ue_4 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_4 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;


% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_12 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_12 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;

% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_24 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_24 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_25 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_25 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_26 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_26 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_27 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_27 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_28 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_28 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_30 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_30 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_31 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_31 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_32 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_32 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_33 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_33 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_34 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_34 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_35 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_35 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_36 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_36 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_37 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_37 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_38 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_38 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_39 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_39 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_40 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_40 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_41 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_41 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_42 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_42 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_43 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_43 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_44 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_44 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_45 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_45 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_46 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_46 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_47 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_47 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_48 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_48 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_49 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_49 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_50 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_50 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_51 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_51 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_52 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_52 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_53 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_53 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_54 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_54 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_55 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_55 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_56 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_56 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_57 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_57 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_58 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_58 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_59 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_59 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_60 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_60 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_61 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_61 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_62 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_62 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_63 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_63 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_64 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_64 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_65 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_65 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_66 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_66 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_67 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_67 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_68 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_68 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_69 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_69 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_70 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_70 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_71 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_71 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_72 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_72 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_73 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_73 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;

% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_74 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_74 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_75 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_75 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_76 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_76 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_77 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_77 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_78 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_78 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_79 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_79 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_80 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_80 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_81 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_81 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_82 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_82 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_83 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_83 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_84 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_84 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_85 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_85 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_86 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_86 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_87 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_87 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_88 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_88 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_89 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_89 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_90 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_90 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_91 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_91 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_92 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_92 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_93 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_93 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_94 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_94 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_95 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_95 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_96 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_96 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_97 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_97 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_98 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_98 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_99 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_99 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_100 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_100 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_101 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_101 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_102 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_102 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_103 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_103 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_104 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;
load Ye_GA_3mix_3_WLS_ue_104 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat;
% Load the nth set of data
load Ye_GA_3mix_3_WLS_ue_105 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + 3;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape([1:3],:,:);
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted([1:3],:,:);
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon([1:3],:,:);
load Ye_GA_3mix_3_WLS_ue_105 recorder_pc_start_repeat
recorder_pc_start_repeat_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat([1:3],:,:);

% Save the final result
clear recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
recorder_qreg_repeat_reshape = recorder_qreg_repeat_reshape_temp; clear recorder_qreg_repeat_reshape_temp
recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted_temp; clear recorder_WLS_start_repeat_sorted_temp
recorder_pc_start_repeat_recon = recorder_pc_start_repeat_recon_temp; clear recorder_pc_start_repeat_recon_temp
recorder_pc_start_repeat = recorder_pc_start_repeat_temp; clear recorder_pc_start_repeat_temp


n_repeatsample = n_repeatsample_temp; clear n_repeatsample_temp
clear j_begin  j_start 
j = n_repeatsample
save Ye_GA_3mix_3_WLS_ue_combined





return;
% The following is temporarily commented 
% The case for piecewise linear estimation for a Laplace EIV
% The variables that we need: 
% taugrid_ue, taugrid, taugrid_midpoint, 
% recorder_qreg_repeat_reshape, recorder_WLS_start_repeat_sorted, 
% recorder_pc_start_repeat_recon,   

% Load tau grid and ncovar
load Ye_GA_L_3_WLS_ue_pl1 taugrid_ue taugrid taugrid_midpoint ncovar

% Load the first set of data
load Ye_GA_L_3_WLS_ue_pl1 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
recorder_qreg_repeat_reshape_temp = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp = recorder_pc_start_repeat_recon;
n_repeatsample_temp = n_repeatsample;

% Load the second set of data
load Ye_GA_L_3_WLS_ue_pl2 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;


% Load the nth set of data
load Ye_GA_L_3_WLS_ue_pl3 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;


% Load the nth set of data
load Ye_GA_L_3_WLS_ue_pl4 recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
j_begin = n_repeatsample_temp + 1;
j_start = n_repeatsample_temp + n_repeatsample;
n_repeatsample_temp = j_start;
recorder_qreg_repeat_reshape_temp([j_begin:j_start],:,:) = recorder_qreg_repeat_reshape;
recorder_WLS_start_repeat_sorted_temp([j_begin:j_start],:,:) = recorder_WLS_start_repeat_sorted;
recorder_pc_start_repeat_recon_temp([j_begin:j_start],:,:) = recorder_pc_start_repeat_recon;


% Save the final result
clear recorder_qreg_repeat_reshape recorder_WLS_start_repeat_sorted recorder_pc_start_repeat_recon n_repeatsample
recorder_qreg_repeat_reshape = recorder_qreg_repeat_reshape_temp; clear recorder_qreg_repeat_reshape_temp
recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted_temp; clear recorder_WLS_start_repeat_sorted_temp
recorder_pc_start_repeat_recon = recorder_pc_start_repeat_recon_temp; clear recorder_pc_start_repeat_recon_temp
n_repeatsample = n_repeatsample_temp; clear n_repeatsample_temp
clear j_begin  j_start 
j = n_repeatsample
save Ye_GA_L_3_WLS_ue_pl_combined











