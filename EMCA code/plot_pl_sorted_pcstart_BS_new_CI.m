% This code plots the full sample piecewise constant and linear results with confidence interval from BS
% Real data
% Written on 10/8/2018
clear; 
close all;

  %%% Part A) Set parameters and load data
    %%% A-1) Full sample result
    % 80 White 10 tau
    % load  pl_Angrist80_sub1_ntau15_white1_result
    % 90 White 10 tau
    % load  pl_Angrist90_sub1_ntau15_white1_result
    % 00 White 15 tau
    % load pl_Angrist00_sub1_ntau15_white1_result 
    % 10 White 15 tau
    load pl_Jacob10_sub1_ntau15_white1_result 
    beta_qreg_full = beta_qreg; clear beta_qreg; 
    beta_WLS_start_sorted_full = beta_WLS_start_sorted; clear beta_WLS_start_sorted
    beta_pc_s_full = beta_pc_s; clear beta_pc_s;
   
    %%% A-2) Load BS data
    % 80 White 10 tau
    % load RData_pl_Angrist80_sub21_ntau15_white100_new
    % 90 White 10 tau
    % load RData_pl_Angrist90_sub21_ntau15_white100_new 
    % 00 White 15 tau
    % load RData_pl_Angrist00_sub21_ntau15_white100_new
    % 10 White 15 tau
    load RData_pl_Jacob10_sub21_ntau15_white20_new
    

  %%% Part B) Preallocation
    beta_qreg_repeat = nan(n_repeatsample,ncovar,ntau);
    beta_WLS_start_sorted_repeat = nan(n_repeatsample,ncovar,ntau);
    beta_pc_s_repeat = nan(n_repeatsample,ncovar,ntau);
 
  
  %%% Part C) Raw estimates (constants + increments) => constructed
  %%% estimates (constants)
   for j_repeat = [1:n_repeatsample]
    % C-1) Load the raw estimates
        % X_mean
        X_mean = X_mean_repeat(:,j_repeat)';
        % qreg 
        recorder_qreg = squeeze(recorder_qreg_repeat(:,:,j_repeat)); 
        % WLS start 
        recorder_WLS_start = recorder_WLS_start_repeat(j_repeat,:);
        % Piecewise constant MLE + Piecewise linear MLE (fmincon)
        recorder_pc_start = recorder_pc_start_repeat(j_repeat,:);
        % Piecewise constant MLE + Piecewise linear MLE (fmincon) + Piecewise linear MLE (GA)
        recorder_fmincon_start = recorder_fmincon_start_repeat(j_repeat,:);
        
    % C-2) Construct the sorted estimates
        % qreg        
        beta_qreg = (reshape(recorder_qreg(1,[1:(ncovar*ntau)]), [ntau, ncovar]))';   
        % WLS start result
        recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau*ncovar),ntau,ncovar);
        [beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,ncovar);
        beta_WLS_start_sorted = beta_WLS_start_sorted';
        % pc start result
        beta_pc_s = reconstruct_beta( (reshape(recorder_pc_start(1,[1:(ncovar*ntau)]), [ntau, ncovar]))'  );
        
      
    % C-3) Save the sorted estimates
        % qreg   
        beta_qreg_repeat(j_repeat,:,:) = beta_qreg;   
        % WLS start result   
        beta_WLS_start_sorted_repeat(j_repeat,:,:)  = beta_WLS_start_sorted;
        % pc start result
        beta_pc_s_repeat(j_repeat,:,:) = beta_pc_s;
        
        clear X_mean beta_qreg beta_WLS_start_sorted beta_pc_s beta_fmincon_s
       
   end
 
   
  %%% Part D) Plotting 
   beta2_pc_s_repreat = squeeze(beta_pc_s_repeat(:,2,:));
    % D-1) Calculate the standard errors
    
    beta_qreg_std = squeeze(std(beta_qreg_repeat));
    beta_WLS_start_sorted_std = squeeze(std(beta_WLS_start_sorted_repeat));
    beta_pc_s_repeat = beta_pc_s_repeat(abs(beta_pc_s_repeat(:,2,end))< 1,:,:); 
  
    beta_pc_s_std = squeeze(std(beta_pc_s_repeat));
    
    beta_pc_s_upper =  prctile(beta2_pc_s_repreat,97.5);
    beta_pc_s_lower =  prctile(beta2_pc_s_repreat,2.5);
  
    
    % D-2) Plot 
    j_covar = 1;
    
    figure; 
    hold on;
    
    output_matrix = zeros(9, 15);
    
    plot(taugrid,beta_qreg_full(j_covar+1,:),'g','LineWidth',1);
    output_matrix(1, :) = beta_qreg_full(j_covar+1,:);
    plot(taugrid,beta_pc_s_full(j_covar+1,:),'r','LineWidth',1);
    output_matrix(3, :) = beta_pc_s_full(j_covar+1,:);
 
    

    
    plot(taugrid,beta_pc_s_upper,'k--');
    plot(taugrid,beta_pc_s_full(j_covar+1,:) + 1.96*beta_pc_s_std(j_covar+1,:),'r--');
    
    
    plot(taugrid,beta_pc_s_lower,'k--');

     

  
    plot(taugrid,beta_pc_s_full(j_covar+1,:) - 1.96*beta_pc_s_std(j_covar+1,:),'r--');
  
   
 
    legend('quantile regression','piecewise linear(fmincon)','CI 2.5%-97.5%','CI 1.96*std'); 
    title(sprintf('Data: %s. # of bootstrap samples: %d',postfixdata, n_repeatsample))
    
    output_matrix = [taugrid;output_matrix]';
    % csvwrite(sprintf('%s.csv',postfixdata),output_matrix,1,0);
    

    
    print('-dpng','-r0',sprintf('BS_pl_CI_%s_%d',postfix,n_repeatsample));
    
    
   
   

   
   
   
    
