% This code plots the full sample piecewise CONSTANT result with confidence interval from BS
% Written on 7/3/2018
clear; 
close all;

  %%% Part A) Set parameters and load data
   %%% A-1) Full sample result
    % 80 White 10 tau
    % load pl_Angrist80_sub1_ntau10_white1_result 
    % 90 White 10 tau
   %  load pl_Angrist90_sub1_ntau10_white1_result 
    % 00 White 15 tau
     load pl_Angrist00_sub1_ntau15_white1_result 
    % load RData_pl_Angrist00_sub1_ntau15_white1
    % 10 White 15 tau
    % load pl_Jacob10_sub1_ntau15_white1_result 
    beta_qreg_full = beta_qreg; clear beta_qreg; 
    beta_WLS_start_sorted_full = beta_WLS_start_sorted; clear beta_WLS_start_sorted
    beta_pc_s_full = beta_pc_s; clear beta_pc_s;
    beta_fmincon_s_full = beta_fmincon_s; clear beta_fmincon_s;
  
    
   %%% A-2) Load BS data
    % 80 White 10 tau
    % load RData_pl_Angrist80_sub21_ntau10_white200c
    % 90 White 10 tau
    % load  RData_pl_Angrist90_sub21_ntau10_white200c
    % 00 White 15 tau
    load RData_pl_Angrist00_sub21_ntau15_white209c
    % 10 White 15 tau
    % load RData_pl_Jacob10_sub21_ntau15_white215c
    
    ntau_qreg = ntau_plus - 2;

    
    
  %%% Part B) Preallocation
    beta_qreg_repeat = nan(n_repeatsample,ncovar,ntau_qreg);
    beta_WLS_start_sorted_repeat = nan(n_repeatsample,ncovar,ntau_qreg);
    beta_pc_s_repeat = nan(n_repeatsample,ncovar,ntau_plus);
    beta_fmincon_s_repeat = nan(n_repeatsample,ncovar,ntau_plus);
  
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
        beta_qreg = (reshape(recorder_qreg(1,[1:(ncovar*ntau_qreg)]), [ntau_qreg, ncovar]))';   
        % WLS start result
        recorder_WLS_start_reshape = reshape(recorder_WLS_start(1:ntau_qreg*ncovar),ntau_qreg,ncovar);
        [beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau_qreg,ncovar);
        beta_WLS_start_sorted = beta_WLS_start_sorted';
        % pc start result
        beta_pc_s = reconstruct_beta( (reshape(recorder_pc_start(1,[1:(ncovar*ntau_plus)]), [ntau_plus, ncovar]))'  );
        % fmincon start result
        beta_fmincon_s = reconstruct_beta( (reshape(recorder_fmincon_start(1,[1:(ncovar*ntau_plus)]), [ntau_plus, ncovar]))'  );
    
    % C-3) Save the sorted estimates
        % qreg   
        beta_qreg_repeat(j_repeat,:,:) = beta_qreg;   
        % WLS start result   
        beta_WLS_start_sorted_repeat(j_repeat,:,:)  = beta_WLS_start_sorted;
        % pc start result
        beta_pc_s_repeat(j_repeat,:,:) = beta_pc_s;
        % fmincon start result
        beta_fmincon_s_repeat(j_repeat,:,:) = beta_fmincon_s;   
        clear X_mean beta_qreg beta_WLS_start_sorted beta_pc_s beta_fmincon_s
       
   end
   
   
  %%% Part D) Plotting 
   beta2_pc_s_repreat = squeeze(beta_pc_s_repeat(:,2,end));
    % D-1) Calculate the standard errors
    
    beta_qreg_std = squeeze(std(beta_qreg_repeat));
    beta_WLS_start_sorted_std = squeeze(std(beta_WLS_start_sorted_repeat));
    beta_pc_s_repeat = beta_pc_s_repeat(abs(beta_pc_s_repeat(:,2,end))< 1,:,:); 
    beta_fmincon_s_repeat = beta_fmincon_s_repeat(abs(beta_fmincon_s_repeat(:,2,end))<1,:,:); 
    beta_pc_s_std = squeeze(std(beta_pc_s_repeat));
    beta_fmincon_s_std = squeeze(std(beta_fmincon_s_repeat));
    
    % D-2) Plot 
    j_covar = 1;
    
    figure; 
    hold on;
    
    plot(taugrid_qreg,beta_qreg_full(j_covar+1,:),'g','LineWidth',1);
    plot(taugrid_qreg,beta_WLS_start_sorted_full(j_covar+1,:),'k','LineWidth',1);
    plot(taugrid_break,beta_pc_s_full(j_covar+1,:),'r','LineWidth',1);
    plot(taugrid_break,beta_fmincon_s_full(j_covar+1,:),'b','LineWidth',1);
    
    
    
    plot(taugrid_qreg,beta_qreg_full(j_covar+1,:) + 1.96*beta_qreg_std(j_covar+1,:),'g--');
    plot(taugrid_qreg,beta_qreg_full(j_covar+1,:) - 1.96*beta_qreg_std(j_covar+1,:),'g--');
    
    plot(taugrid_qreg,beta_WLS_start_sorted_full(j_covar+1,:) + 1.96*beta_WLS_start_sorted_std(j_covar+1,:),'k--');
    plot(taugrid_qreg,beta_WLS_start_sorted_full(j_covar+1,:) - 1.96*beta_WLS_start_sorted_std(j_covar+1,:),'k--');
    
    plot(taugrid_break,beta_pc_s_full(j_covar+1,:) + 1.96*beta_pc_s_std(j_covar+1,:),'r--');
    plot(taugrid_break,beta_pc_s_full(j_covar+1,:) - 1.96*beta_pc_s_std(j_covar+1,:),'r--');
    
    plot(taugrid_break,beta_fmincon_s_full(j_covar+1,:) + 1.96*beta_fmincon_s_std(j_covar+1,:),'b--');
    plot(taugrid_break,beta_fmincon_s_full(j_covar+1,:) - 1.96*beta_fmincon_s_std(j_covar+1,:),'b--');   
    
    legend('quantile regression','piecewise constant','piecewise linear(fmincon)', 'piecewise linear(GA)'); 
    title(sprintf('Data: %s. # of bootstrap samples: %d',postfixdata, n_repeatsample))
    print('-dpng','-r0',sprintf('BS_pl_%s_%d',postfix,n_repeatsample));
    
    
   
   

   
   
   
    
