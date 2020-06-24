% This code plots the piecewise LINEAR real data result including the two
% end points
% The code is based on plot_pl_sorted_GA_pcstart.m for the results without
% the end points
% Updated on 7/28/2018
% Haoyang Liu

clear
close all 

  %%% Part A) Set Parameters and Load Data
    % 80 White 10 tau
    load  RData_pl_end_Angrist80_sub1_ntau9_white1
    % 90 White 10 tau
    % load RData_pl_end_Angrist90_sub1_ntau10_white1 
    % 00 White 15 tau
    % load RData_pl_end_Angrist00_sub1_ntau15_white1 
    % 10 White 15 tau
    % load  RData_pl_end_Jacob10_sub1_ntau14_white1
    % The following adjustments also need to be made for white subsamples
    taugrid = taugrid_pl; 
    % Here ntau is the number of segments, instead of the number of barriers
    ntau = ntau + 1;
    postfixdata = 'Angrist 80 (White)'; % Jacob 10
  %%% End of Part A)

  %%% Part B) Reconstruct beta
    % pc start result
    beta_pc_s = reconstruct_beta( (reshape(recorder_pc_start_repeat(1,[1:(ncovar*ntau_plus)]), [ntau_plus, ncovar]))'  );

    beta_sorted_pc_s = nan(ncovar, 2*ntau);
  
  %%% End of Part B)

  %%% Part C) Create "taugrid_break_sorted", which repeats each middle entry twice 
    step_original = taugrid(3) - taugrid(2);
    ntau_qreg = ntau_plus - 2;
    % first entry is 0 instead of 0+0.00001
    taugrid_break = [0:step_original:1];
    taugrid_break_sorted = repeat_HL(taugrid_break,2);
    % Repeating it twice gives us end points
    taugrid_break_sorted = taugrid_break_sorted(1:end-1);
  %%% End of Part C)

  %%% Part D) Fill in the parametes for qreg
    beta_qreg = (reshape(recorder_qreg(1,[1:(ncovar*ntau_plus)]), [ntau_plus, ncovar]))';
  %%% End of Part D)

  %%% Part E) Sorting 
    beta_mid_pc_s = (beta_pc_s(:, [1:end-1]) + beta_pc_s(:, [2:end]))/2;
    fval_tau_pc_s = X_mean_repeat'*beta_mid_pc_s;
    [V_pc_s,I_pc_s] = sort(fval_tau_pc_s);
 
    beta_sorted_pc_s(:,[1:2:end]) = beta_pc_s(:,I_pc_s);
    beta_sorted_pc_s(:,[2:2:end]) = beta_pc_s(:,I_pc_s+1);
  %%% End of Part E)


  %%% Part F) Plotting
    for j_covar = [1]
    
        figure; hold on;
    
        plot(taugrid_qreg,beta_qreg(j_covar+1,:),'g');
        plot(taugrid_qreg,beta_WLS_start_sorted(j_covar+1,:),'k');
        plot(taugrid_break,beta_pc_s(j_covar+1,:),'r');
        
    
        legend('quantile regression','piecewise constant','piecewise linear');
        title(sprintf('Data: %s. Beta %d',postfixdata, j_covar))
        print('-dpng','-r0',sprintf('beta%d_%s_end',j_covar,postfix));
    end
    
