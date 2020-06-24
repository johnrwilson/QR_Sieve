% This code combines all bootstrap results of piecewise linear for 2000 and 2010
% 7/26/2018
clear 
close all 

%%% 1980
 % Load the first dataset
    load RData_pl_Angrist80_sub21_ntau10_white50 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = X_mean_repeat';
    recorder_qreg_repeat_all = squeeze(recorder_qreg_repeat)';
    recorder_WLS_start_repeat_all = recorder_WLS_start_repeat;
    recorder_pc_start_repeat_all = recorder_pc_start_repeat;
    recorder_fmincon_start_repeat_all = recorder_fmincon_start_repeat;
    n_repeatsample_all = j_repeatsample;
    clear j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat

  % Load the second dataset
    load RData_pl_Angrist80_sub21_ntau10_white150 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = [X_mean_repeat_all; X_mean_repeat(:, [1:j_repeatsample])'];
    recorder_qreg_repeat_all = [recorder_qreg_repeat_all; squeeze(recorder_qreg_repeat(:,:,[1:j_repeatsample]))'];
    recorder_WLS_start_repeat_all = [recorder_WLS_start_repeat_all; recorder_WLS_start_repeat([1:j_repeatsample],:)];
    recorder_pc_start_repeat_all = [recorder_pc_start_repeat_all; recorder_pc_start_repeat([1:j_repeatsample],:)];
    recorder_fmincon_start_repeat_all = [recorder_fmincon_start_repeat_all; recorder_fmincon_start_repeat([1:j_repeatsample],:)];
    n_repeatsample_all = n_repeatsample_all + j_repeatsample;
    clear j_repeatsample X_mean_repeat  recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
  
  % Put the combined datasets into the right format
    X_mean_repeat = X_mean_repeat_all'; clear X_mean_repeat_all
    recorder_qreg_repeat(1,:,:) = recorder_qreg_repeat_all'; clear recorder_qreg_repeat_all
    recorder_WLS_start_repeat = recorder_WLS_start_repeat_all; clear recorder_WLS_start_repeat_all
    recorder_pc_start_repeat = recorder_pc_start_repeat_all; clear recorder_pc_start_repeat_all
    recorder_fmincon_start_repeat = recorder_fmincon_start_repeat_all; clear recorder_fmincon_start_repeat_all
    n_repeatsample = n_repeatsample_all; clear n_repeatsample_all

  % Load other constants
    load RData_pl_Angrist80_sub21_ntau10_white50 ntau_plus ncovar taugrid_pl

    save RData_pl_Angrist80_sub21_ntau10_white200c

%%% 1990 
  % Load the first dataset
    load RData_pl_Angrist90_sub21_ntau10_white50 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = X_mean_repeat';
    recorder_qreg_repeat_all = squeeze(recorder_qreg_repeat)';
    recorder_WLS_start_repeat_all = recorder_WLS_start_repeat;
    recorder_pc_start_repeat_all = recorder_pc_start_repeat;
    recorder_fmincon_start_repeat_all = recorder_fmincon_start_repeat;
    n_repeatsample_all = j_repeatsample;
    clear j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat

  % Load the second dataset
    load RData_pl_Angrist90_sub21_ntau10_white150 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = [X_mean_repeat_all; X_mean_repeat(:, [1:j_repeatsample])'];
    recorder_qreg_repeat_all = [recorder_qreg_repeat_all; squeeze(recorder_qreg_repeat(:,:,[1:j_repeatsample]))'];
    recorder_WLS_start_repeat_all = [recorder_WLS_start_repeat_all; recorder_WLS_start_repeat([1:j_repeatsample],:)];
    recorder_pc_start_repeat_all = [recorder_pc_start_repeat_all; recorder_pc_start_repeat([1:j_repeatsample],:)];
    recorder_fmincon_start_repeat_all = [recorder_fmincon_start_repeat_all; recorder_fmincon_start_repeat([1:j_repeatsample],:)];
    n_repeatsample_all = n_repeatsample_all + j_repeatsample;
    clear j_repeatsample X_mean_repeat  recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
  
  % Put the combined datasets into the right format
    X_mean_repeat = X_mean_repeat_all'; clear X_mean_repeat_all
    recorder_qreg_repeat(1,:,:) = recorder_qreg_repeat_all'; clear recorder_qreg_repeat_all
    recorder_WLS_start_repeat = recorder_WLS_start_repeat_all; clear recorder_WLS_start_repeat_all
    recorder_pc_start_repeat = recorder_pc_start_repeat_all; clear recorder_pc_start_repeat_all
    recorder_fmincon_start_repeat = recorder_fmincon_start_repeat_all; clear recorder_fmincon_start_repeat_all
    n_repeatsample = n_repeatsample_all; clear n_repeatsample_all

  % Load other constants
    load RData_pl_Angrist90_sub21_ntau10_white50 ntau_plus ncovar taugrid_pl

    save RData_pl_Angrist90_sub21_ntau10_white200c
    
    
    

%%% 2000 
  % Load the first dataset
    load RData_pl_Angrist00_sub21_ntau15_white20 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = X_mean_repeat';
    recorder_qreg_repeat_all = squeeze(recorder_qreg_repeat)';
    recorder_WLS_start_repeat_all = recorder_WLS_start_repeat;
    recorder_pc_start_repeat_all = recorder_pc_start_repeat;
    recorder_fmincon_start_repeat_all = recorder_fmincon_start_repeat;
    n_repeatsample_all = j_repeatsample;
    clear j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat

  % Load the second dataset
    load RData_pl_Angrist00_sub21_ntau15_white35 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = [X_mean_repeat_all; X_mean_repeat(:, [1:j_repeatsample])'];
    recorder_qreg_repeat_all = [recorder_qreg_repeat_all; squeeze(recorder_qreg_repeat(:,:,[1:j_repeatsample]))'];
    recorder_WLS_start_repeat_all = [recorder_WLS_start_repeat_all; recorder_WLS_start_repeat([1:j_repeatsample],:)];
    recorder_pc_start_repeat_all = [recorder_pc_start_repeat_all; recorder_pc_start_repeat([1:j_repeatsample],:)];
    recorder_fmincon_start_repeat_all = [recorder_fmincon_start_repeat_all; recorder_fmincon_start_repeat([1:j_repeatsample],:)];
    n_repeatsample_all = n_repeatsample_all + j_repeatsample;
    clear j_repeatsample X_mean_repeat  recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat

  % Load the third one
    load RData_pl_Angrist00_sub21_ntau15_white44 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = [X_mean_repeat_all; X_mean_repeat(:, [1:j_repeatsample])'];
    recorder_qreg_repeat_all = [recorder_qreg_repeat_all; squeeze(recorder_qreg_repeat(:,:,[1:j_repeatsample]))'];
    recorder_WLS_start_repeat_all = [recorder_WLS_start_repeat_all; recorder_WLS_start_repeat([1:j_repeatsample],:)];
    recorder_pc_start_repeat_all = [recorder_pc_start_repeat_all; recorder_pc_start_repeat([1:j_repeatsample],:)];
    recorder_fmincon_start_repeat_all = [recorder_fmincon_start_repeat_all; recorder_fmincon_start_repeat([1:j_repeatsample],:)];
    n_repeatsample_all = n_repeatsample_all + j_repeatsample;
    clear j_repeatsample X_mean_repeat  recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
  
    
  % Load the fourth one
    load RData_pl_Angrist00_sub21_ntau15_white110 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = [X_mean_repeat_all; X_mean_repeat(:, [1:j_repeatsample])'];
    recorder_qreg_repeat_all = [recorder_qreg_repeat_all; squeeze(recorder_qreg_repeat(:,:,[1:j_repeatsample]))'];
    recorder_WLS_start_repeat_all = [recorder_WLS_start_repeat_all; recorder_WLS_start_repeat([1:j_repeatsample],:)];
    recorder_pc_start_repeat_all = [recorder_pc_start_repeat_all; recorder_pc_start_repeat([1:j_repeatsample],:)];
    recorder_fmincon_start_repeat_all = [recorder_fmincon_start_repeat_all; recorder_fmincon_start_repeat([1:j_repeatsample],:)];
    n_repeatsample_all = n_repeatsample_all + j_repeatsample;
    clear j_repeatsample X_mean_repeat  recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat

    
  % Put the combined datasets into the right format
    X_mean_repeat = X_mean_repeat_all'; clear X_mean_repeat_all
    recorder_qreg_repeat(1,:,:) = recorder_qreg_repeat_all'; clear recorder_qreg_repeat_all
    recorder_WLS_start_repeat = recorder_WLS_start_repeat_all; clear recorder_WLS_start_repeat_all
    recorder_pc_start_repeat = recorder_pc_start_repeat_all; clear recorder_pc_start_repeat_all
    recorder_fmincon_start_repeat = recorder_fmincon_start_repeat_all; clear recorder_fmincon_start_repeat_all
    n_repeatsample = n_repeatsample_all; clear n_repeatsample_all

  % Load other constants
    load RData_pl_Angrist00_sub21_ntau15_white20 ntau_plus ncovar taugrid_pl

    save RData_pl_Angrist00_sub21_ntau15_white209c


%%% 2010 
  % Load the first dataset
    load RData_pl_Jacob10_sub21_ntau15_white20 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = X_mean_repeat';
    recorder_qreg_repeat_all = squeeze(recorder_qreg_repeat)';
    recorder_WLS_start_repeat_all = recorder_WLS_start_repeat;
    recorder_pc_start_repeat_all = recorder_pc_start_repeat;
    recorder_fmincon_start_repeat_all = recorder_fmincon_start_repeat;
    n_repeatsample_all = j_repeatsample;
    clear j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat

  % Load the second dataset
    load RData_pl_Jacob10_sub21_ntau15_white38 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = [X_mean_repeat_all; X_mean_repeat(:, [1:j_repeatsample])'];
    recorder_qreg_repeat_all = [recorder_qreg_repeat_all; squeeze(recorder_qreg_repeat(:,:,[1:j_repeatsample]))'];
    recorder_WLS_start_repeat_all = [recorder_WLS_start_repeat_all; recorder_WLS_start_repeat([1:j_repeatsample],:)];
    recorder_pc_start_repeat_all = [recorder_pc_start_repeat_all; recorder_pc_start_repeat([1:j_repeatsample],:)];
    recorder_fmincon_start_repeat_all = [recorder_fmincon_start_repeat_all; recorder_fmincon_start_repeat([1:j_repeatsample],:)];
    n_repeatsample_all = n_repeatsample_all + j_repeatsample;
    clear j_repeatsample X_mean_repeat  recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat

  % Load the third one
    load RData_pl_Jacob10_sub21_ntau15_white47 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = [X_mean_repeat_all; X_mean_repeat(:, [1:j_repeatsample])'];
    recorder_qreg_repeat_all = [recorder_qreg_repeat_all; squeeze(recorder_qreg_repeat(:,:,[1:j_repeatsample]))'];
    recorder_WLS_start_repeat_all = [recorder_WLS_start_repeat_all; recorder_WLS_start_repeat([1:j_repeatsample],:)];
    recorder_pc_start_repeat_all = [recorder_pc_start_repeat_all; recorder_pc_start_repeat([1:j_repeatsample],:)];
    recorder_fmincon_start_repeat_all = [recorder_fmincon_start_repeat_all; recorder_fmincon_start_repeat([1:j_repeatsample],:)];
    n_repeatsample_all = n_repeatsample_all + j_repeatsample;
    clear j_repeatsample X_mean_repeat  recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
  
  % Load the fourth one
    load RData_pl_Jacob10_sub21_ntau15_white110 j_repeatsample X_mean_repeat recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    X_mean_repeat_all = [X_mean_repeat_all; X_mean_repeat(:, [1:j_repeatsample])'];
    recorder_qreg_repeat_all = [recorder_qreg_repeat_all; squeeze(recorder_qreg_repeat(:,:,[1:j_repeatsample]))'];
    recorder_WLS_start_repeat_all = [recorder_WLS_start_repeat_all; recorder_WLS_start_repeat([1:j_repeatsample],:)];
    recorder_pc_start_repeat_all = [recorder_pc_start_repeat_all; recorder_pc_start_repeat([1:j_repeatsample],:)];
    recorder_fmincon_start_repeat_all = [recorder_fmincon_start_repeat_all; recorder_fmincon_start_repeat([1:j_repeatsample],:)];
    n_repeatsample_all = n_repeatsample_all + j_repeatsample;
    clear j_repeatsample X_mean_repeat  recorder_qreg_repeat recorder_WLS_start_repeat recorder_pc_start_repeat recorder_fmincon_start_repeat
    
    
  % Put the combined datasets into the right format
    X_mean_repeat = X_mean_repeat_all'; clear X_mean_repeat_all
    recorder_qreg_repeat(1,:,:) = recorder_qreg_repeat_all'; clear recorder_qreg_repeat_all
    recorder_WLS_start_repeat = recorder_WLS_start_repeat_all; clear recorder_WLS_start_repeat_all
    recorder_pc_start_repeat = recorder_pc_start_repeat_all; clear recorder_pc_start_repeat_all
    recorder_fmincon_start_repeat = recorder_fmincon_start_repeat_all; clear recorder_fmincon_start_repeat_all
    n_repeatsample = n_repeatsample_all; clear n_repeatsample_all

  % Load other constants
    load RData_pl_Jacob10_sub21_ntau15_white20 ntau_plus ncovar taugrid_pl

    % save RData_pl_Jacob10_sub21_ntau15_white105c
    save RData_pl_Jacob10_sub21_ntau15_white215c







 
 
 
