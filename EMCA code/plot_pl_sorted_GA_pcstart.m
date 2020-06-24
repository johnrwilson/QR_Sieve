% This code plots the piecewise LINEAR full sample real data result
% What makes piecewise LINEAR different from piecewise CONSTANT are:
% 1) In sorting, we need to move around the segments, instead of individual
% points (in part C) and part E)) 
% 2) The parameters for piecewise linear are constants plus increments,
% instead of all constants (in part B))
% It can plot the old full sample results, the new subsample, and the new white subsample results 
% This is the new version that allows qreg grid to be shorter than fmincon grid
% It can also plot sorted estimates (currently commented out) 
% Updated on 7/3/2018
% Haoyang Liu

clear
close all 

  %%% Part A) Set Parameters and Load Data
    %postfix = 'Angrist80' 'Jacob10'; (This should be for old full sample result. Didn't check if this still works on 7/3)
    postfix = 'Jacob10'
    postfixdata = 'Jacob 10' % Jacob 10

    
    % 80 White 10 tau
    % load RData_pl_Angrist80_sub1_ntau15_white1 
    % 90 White 10 tau
    % load RData_pl_Angrist90_sub1_ntau19_white1 
    % 00 White 15 tau
    % load RData_pl_Angrist00_sub1_ntau17_white1 
    % 10 White 15 tau
    load RData_pl_Jacob10_sub1_ntau15_white1 
    % The following adjustments also need to be made for white subsamples
    taugrid = taugrid_pl; 
    % Here ntau is the number of segments, instead of the number of barriers
    ntau = ntau + 1;
  
    
    % ntau_s = '34_pcstart_s'; (Didn't check if the following still works on 7/3)
    % 80
    % load (sprintf('GA_%s_sub65022_ntau%s',postfix,ntau_s))
    % 90
    % load (sprintf('GA_%s_sub86784_ntau%s',postfix,ntau_s))
    % 00 
    % load (sprintf('GA_%s_sub97396_ntau%s',postfix,ntau_s))
    % 10
    % load (sprintf('GA_%s_sub106624_ntau%s',postfix,ntau_s))
  %%% End of Part A)

  %%% Part B) Reconstruct beta
    % pc start result
    beta_pc_s = reconstruct_beta( (reshape(recorder_pc_start_repeat(1,[1:(ncovar*ntau_plus)]), [ntau_plus, ncovar]))'  );
    % fmincon start result
    beta_fmincon_s = reconstruct_beta( (reshape(recorder_fmincon_start_repeat(1,[1:(ncovar*ntau_plus)]), [ntau_plus, ncovar]))'  );

    beta_sorted_pc_s = nan(ncovar, 2*ntau);
    beta_sorted_fmincon_s = nan(ncovar, 2*ntau);
  %%% End of Part B)

  %%% Part C) Create "taugrid_break_sorted", which repeats each middle entry twice 
    step_original = taugrid(3) - taugrid(2);
    taugrid_qreg = taugrid([2:end-1]);
    ntau_qreg = ntau_plus - 2;
    % first entry is 0 instead of 0+0.00001
    taugrid_break = [0:step_original:1];
    taugrid_break_sorted = repeat_HL(taugrid_break,2);
    % Repeating it twice gives us end points
    taugrid_break_sorted = taugrid_break_sorted(1:end-1);
  %%% End of Part C)

  %%% Part D) Fill in the parametes for qreg
    beta_qreg = (reshape(recorder_qreg(1,[1:(ncovar*ntau_qreg)]), [ntau_qreg, ncovar]))';
  %%% End of Part D)

  %%% Part E) Sorting 
    beta_mid_pc_s = (beta_pc_s(:, [1:end-1]) + beta_pc_s(:, [2:end]))/2;
    beta_mid_fmincon_s = (beta_fmincon_s(:, [1:end-1]) + beta_fmincon_s(:, [2:end]))/2;

    fval_tau_pc_s = X_mean_repeat'*beta_mid_pc_s;
    fval_tau_fmincon_s = X_mean_repeat'*beta_mid_fmincon_s;

    [V_pc_s,I_pc_s] = sort(fval_tau_pc_s);
    [V_fmincon_s,I_fmincon_s] = sort(fval_tau_fmincon_s);

    beta_sorted_pc_s(:,[1:2:end]) = beta_pc_s(:,I_pc_s);
    beta_sorted_pc_s(:,[2:2:end]) = beta_pc_s(:,I_pc_s+1);
    beta_sorted_fmincon_s(:,[1:2:end]) = beta_fmincon_s(:,I_fmincon_s);
    beta_sorted_fmincon_s(:,[2:2:end]) = beta_fmincon_s(:,I_fmincon_s+1);
  %%% End of Part E)


  %%% Part F) Plotting
    for j_covar = [1]
    
        figure; hold on;
    
        plot(taugrid_qreg,beta_qreg(j_covar+1,:),'g');
        plot(taugrid_qreg,beta_WLS_start_sorted(j_covar+1,:),'k');
        plot(taugrid_break,beta_pc_s(j_covar+1,:),'r');
        % plot(taugrid_break_sorted,beta_sorted_pc_s(j_covar+1,:),'r--');
        
        % plot(taugrid_break,beta_fmincon_s(j_covar+1,:),'b');
        % plot(taugrid_break_sorted,beta_sorted_fmincon_s(j_covar+1,:),'b--');
    
        legend('quantile regression','piecewise constant','piecewise linear(fmincon)');% , 'piecewise linear(GA)'
        title(sprintf('Data: %s. Number of knots %d',postfixdata, (ntau-1)))
        print('-dpng','-r0',sprintf('beta%d_%s',(ntau-1),postfix));
    end
    
    % save pl_Angrist00_sub1_ntau15_white1_result taugrid_qreg taugrid_break beta_qreg beta_WLS_start_sorted beta_pc_s beta_fmincon_s postfixdata postfix


