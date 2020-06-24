% This file combines simualted data for new

% 11/14/2018
clear


%%% T
% Load the first dataset
    load Ye_GA_T_3_par_start_3 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = exit_recorder_pc_start;
    x1matrix_all = x1matrix;
    x2matrix_all = x2matrix;
    recorder_qreg_all = recorder_qreg;
    recorder_qreg_start_all = recorder_qreg_start;
    recorder_pc_start_all = recorder_pc_start;
    clear exit_recorder_pc_start  x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
   

    temp = [1:length(exit_recorder_pc_start_all)];
    
    temp = temp.*(exit_recorder_pc_start_all~=0);
    temp = temp(temp>0);
    exit_recorder_pc_start_all = exit_recorder_pc_start_all(temp);
    x1matrix_all = x1matrix_all(temp,:);
    x2matrix_all = x2matrix_all(temp,:);
    recorder_qreg_all = recorder_qreg_all(temp,:);
    recorder_pc_start_all = recorder_pc_start_all(temp,:);
    recorder_qreg_start_all = recorder_qreg_start_all(temp,:);
    

  % Put the combined datasets into the right format
    exit_recorder_pc_start = exit_recorder_pc_start_all([1:100]); clear exit_recorder_pc_start_all
    x1matrix = x1matrix_all([1:100],:); clear x1matrix_all
    x2matrix = x2matrix_all([1:100],:); clear x2matrix_all
    recorder_qreg = recorder_qreg_all([1:100],:); clear recorder_qreg_all
    recorder_qreg_start = recorder_qreg_start_all([1:100],:); clear recorder_qreg_start_all
    recorder_pc_start = recorder_pc_start_all([1:100],:); clear recorder_pc_start_all

  % Load other constants
    load Ye_GA_T_3_par_start_3 ntau taugrid beta0_true beta1_true beta2_true

    save Ye_GA_T_3_par_start_combined_1
    
    
    
return;

%%% L
% Load the first dataset
    load Ye_GA_L_3_par_start_new exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = exit_recorder_pc_start;
    x1matrix_all = x1matrix;
    x2matrix_all = x2matrix;
    recorder_qreg_all = recorder_qreg;
    recorder_qreg_start_all = recorder_qreg_start;
    recorder_pc_start_all = recorder_pc_start;
    clear exit_recorder_pc_start  x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    
    
% Load the second dataset
    load Ye_GA_L_3_par_start_1 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    
% Load the third dataset
    load Ye_GA_L_3_par_start_2 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start


 % Load the fourth dataset
    load Ye_GA_L_3_par_start_3 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start

  
 % Load the fifth dataset
    load Ye_GA_L_3_par_start_4 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start

  
 % Load the sixth dataset
    load Ye_GA_L_3_par_start_5 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start

  % Put the combined datasets into the right format
    exit_recorder_pc_start = exit_recorder_pc_start_all([1:100]); clear exit_recorder_pc_start_all
    x1matrix = x1matrix_all([1:100],:); clear x1matrix_all
    x2matrix = x2matrix_all([1:100],:); clear x2matrix_all
    recorder_qreg = recorder_qreg_all([1:100],:); clear recorder_qreg_all
    recorder_qreg_start = recorder_qreg_start_all([1:100],:); clear recorder_qreg_start_all
    recorder_pc_start = recorder_pc_start_all([1:100],:); clear recorder_pc_start_all

  % Load other constants
    load Ye_GA_L_3_par_start_new ntau taugrid beta0_true beta1_true beta2_true

    save Ye_GA_L_3_par_start_combined
    
    
%%% mixture of three normal 
% Load the first dataset
    load Ye_GA_3mix_3_par_start_3 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = exit_recorder_pc_start;
    x1matrix_all = x1matrix;
    x2matrix_all = x2matrix;
    recorder_qreg_all = recorder_qreg;
    recorder_qreg_start_all = recorder_qreg_start;
    recorder_pc_start_all = recorder_pc_start;
    clear exit_recorder_pc_start  x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    
    
% Load the second dataset
    load Ye_GA_3mix_3_par_start_3_san exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    
% Load the third dataset
    load Ye_GA_3mix_3_par_start_3_2 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start


 % Load the fourth dataset
    load Ye_GA_3mix_3_par_start_3_3 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start

  
 % Load the fifth dataset
    load Ye_GA_3mix_3_par_start_3_4 exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start
    exit_recorder_pc_start_all = [exit_recorder_pc_start_all,exit_recorder_pc_start];
    x1matrix_all = [x1matrix_all; x1matrix];
    x2matrix_all = [x2matrix_all; x2matrix];
    recorder_qreg_all = [recorder_qreg_all; recorder_qreg];
    recorder_qreg_start_all = [recorder_qreg_start_all; recorder_qreg_start];
    recorder_pc_start_all = [recorder_pc_start_all; recorder_pc_start];
    clear exit_recorder_pc_start x1matrix x2matrix recorder_qreg recorder_qreg_start recorder_pc_start

  

  % Put the combined datasets into the right format
    exit_recorder_pc_start = exit_recorder_pc_start_all([1:100]); clear exit_recorder_pc_start_all
    x1matrix = x1matrix_all([1:100],:); clear x1matrix_all
    x2matrix = x2matrix_all([1:100],:); clear x2matrix_all
    recorder_qreg = recorder_qreg_all([1:100],:); clear recorder_qreg_all
    recorder_qreg_start = recorder_qreg_start_all([1:100],:); clear recorder_qreg_start_all
    recorder_pc_start = recorder_pc_start_all([1:100],:); clear recorder_pc_start_all

  % Load other constants
    load Ye_GA_L_3_par_start_new ntau taugrid beta0_true beta1_true beta2_true

    save Ye_GA_3mix_3_par_start_combined   
    
    
    
    
    
    
    
