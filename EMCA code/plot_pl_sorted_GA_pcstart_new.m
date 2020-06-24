% This code plots the piecewise LINEAR full sample real data result
% It is the new version that accommodate a shorter tau grid than before
% What makes piecewise LINEAR different from piecewise CONSTANT are:
% 1) In sorting, we need to move around the segments, instead of individual
% points (in part C) and part E)) 
% 2) The parameters for piecewise linear are constants plus increments,
% instead of all constants (in part B))
% It can plot the old full sample results, the new subsample, and the new white subsample results 
% This is the new version that allows qreg grid to be shorter than fmincon grid
% It can also plot sorted estimates (currently commented out) 
% Updated on 8/17/2018
% Haoyang Liu

clear
close all 

  %%% Part A) Set Parameters and Load Data
    %postfix = 'Angrist80' 'Jacob10'; (This should be for old full sample result. Didn't check if this still works on 7/3)
    postfix = 'Angrist00'
    postfixdata = 'Angrist 00' %   'Angrist 90' % 
    
    % 80 White 15 tau
    % load RData_pl_Angrist80_sub1_ntau15_white1_new 
    % 90 White 15 tau
    % load RData_pl_Angrist90_sub1_ntau15_white1_new
    % 00 White 15 tau
    load RData_pl_Angrist00_sub1_ntau15_white1_new 
    % 10 White 15 tau
   % load RData_pl_Jacob10_sub1_ntau15_white1_new    

  %%% End of Part A)

  %%% Part B) Reconstruct beta
    % pc start result
    beta_pc_s = reconstruct_beta( (reshape(recorder_pc_start_repeat(1,[1:(ncovar*ntau)]), [ntau, ncovar]))'  );
    % fmincon start result (Commented because the current version does not use 
    % beta_fmincon_s = reconstruct_beta( (reshape(recorder_fmincon_start_repeat(1,[1:(ncovar*ntau)]), [ntau, ncovar]))'  );
    beta_sorted_pc_s = nan(ncovar, 2*(ntau-1));
    % beta_sorted_fmincon_s = nan(ncovar, 2*(ntau-1));
  %%% End of Part B)

  %%% Part C) Create "taugrid_break_sorted", which repeats each middle entry twice 
    taugrid_break_sorted = repeat_HL(taugrid,2);
    % Repeating it twice gives us end points
    taugrid_break_sorted = taugrid_break_sorted(1:end-1);
  %%% End of Part C)

  %%% Part D) Fill in the parametes for qreg
    beta_qreg = (reshape(recorder_qreg(1,[1:(ncovar*ntau)]), [ntau, ncovar]))';
  %%% End of Part D)

  %%% Part E) Sorting 
    beta_mid_pc_s = (beta_pc_s(:, [1:end-1]) + beta_pc_s(:, [2:end]))/2;
    % beta_mid_fmincon_s = (beta_fmincon_s(:, [1:end-1]) + beta_fmincon_s(:, [2:end]))/2;

    fval_tau_pc_s = X_mean_repeat'*beta_mid_pc_s;
    % fval_tau_fmincon_s = X_mean_repeat'*beta_mid_fmincon_s;

    [V_pc_s,I_pc_s] = sort(fval_tau_pc_s);
    % [V_fmincon_s,I_fmincon_s] = sort(fval_tau_fmincon_s);

    beta_sorted_pc_s(:,[1:2:end]) = beta_pc_s(:,I_pc_s);
    beta_sorted_pc_s(:,[2:2:end]) = beta_pc_s(:,I_pc_s+1);
    % beta_sorted_fmincon_s(:,[1:2:end]) = beta_fmincon_s(:,I_fmincon_s);
    % beta_sorted_fmincon_s(:,[2:2:end]) = beta_fmincon_s(:,I_fmincon_s+1);
  %%% End of Part E)


  %%% Part F) Plotting
    for j_covar = [1]
    
        figure; hold on;
    
        plot(taugrid,beta_qreg(j_covar+1,:),'g');
        plot(taugrid,beta_WLS_start_sorted(j_covar+1,:),'k');
        plot(taugrid,beta_pc_s(j_covar+1,:),'r');
        % plot(taugrid_break_sorted,beta_sorted_pc_s(j_covar+1,:),'r--');
        
        % plot(taugrid_break,beta_fmincon_s(j_covar+1,:),'b');
        % plot(taugrid_break_sorted,beta_sorted_fmincon_s(j_covar+1,:),'b--');
    
        legend('quantile regression','piecewise constant','piecewise linear(fmincon)','Location','northeast');% , 'piecewise linear(GA)'
        title(sprintf('Data: %s. Number of knots %d',postfixdata, (ntau)))
        print('-dpng','-r0',sprintf('beta%d_%s_new',(ntau),postfix));
    end
    
    save pl_Angrist80_sub1_ntau15_white1_result taugrid  beta_qreg beta_WLS_start_sorted beta_pc_s postfixdata postfix


