% This code plots GA real data results
% This is the new version that allows qreg grid to be shorter than fmincon
% grid
% It plots sorted estimates
% It doesn't let sorting done at a finer grid 
% 4/13/2018
% Haoyang Liu

clear
close all 

%%% Part A) Set Parameters and Load Data
postfix = 'Angrist80';
error = 'Mixture of 3';
% ntau_s = '19_fminconstart'; (This option is no longer used). 
% ntau_s = '20_fminconstart_5';
% ntau_s = '20_fminconstart_white_5';
ntau_s = '20_fminconstart_black_5';
% 80
% load (sprintf('GA_%s_sub65022_ntau%s',postfix,ntau_s))
%90
% load (sprintf('GA_%s_sub86784_ntau%s',postfix,ntau_s))
%00 
% load (sprintf('GA_%s_sub97396_ntau%s',postfix,ntau_s))
%10
% load (sprintf('GA_%s_sub106624_ntau%s',postfix,ntau_s))
% 80 White
% load (sprintf('GA_%s_sub60050_ntau%s',postfix,ntau_s))
% 90 White
% load (sprintf('GA_%s_sub80114_ntau%s',postfix,ntau_s))
% 00 White
% load (sprintf('GA_%s_sub90200_ntau%s',postfix,ntau_s))
% 10 White
% load (sprintf('GA_%s_sub90200_ntau%s',postfix,ntau_s))
% 80 Black
load (sprintf('GA_%s_sub4972_ntau%s',postfix,ntau_s))
% 90 Black
% load (sprintf('GA_%s_sub6670_ntau%s',postfix,ntau_s))
% 00 Black
% load (sprintf('GA_%s_sub7196_ntau%s',postfix,ntau_s))
% 10 Black
% load (sprintf('GA_%s_sub8310_ntau%s',postfix,ntau_s))
postfixdata = 'Angrist 80 (Black)'
%%% End of Part A)

%%% Part B) Reconstruct beta
% qreg start result
beta_qreg_s = reconstruct_beta( (reshape(recorder_qreg_start_repeat(1,[1:(ncovar*ntau_plus)]), [ntau_plus, ncovar]))'  );
% fmincon start result
beta_fmincon_s = reconstruct_beta( (reshape(recorder_fmincon_start_repeat(1,[1:(ncovar*ntau_plus)]), [ntau_plus, ncovar]))'  );

beta_sorted_qreg_s = nan(ncovar, 2*ntau);
beta_sorted_fmincon_s = nan(ncovar, 2*ntau);

%%% Part C) Create "taugrid_break_sorted", which repeats each middle entry twice 
step_original = taugrid_qreg(3) - taugrid_qreg(2);
% first entry is 0 instead of 0+0.00001
taugrid_break = [0:step_original:1];
taugrid_break_sorted = repeat_HL(taugrid_break,2);
% Repeating it twice gives us end points
taugrid_break_sorted = taugrid_break_sorted(1:end-1);

%%% Part D) Fill in the parametes for qreg
beta_qreg = (reshape(recorder_qreg(1,[1:(ncovar*ntau_qreg)]), [ntau_qreg, ncovar]))';

%%% Part E) Sorting 
beta_mid_qreg_s = (beta_qreg_s(:, [1:end-1]) + beta_qreg_s(:, [2:end]))/2;
beta_mid_fmincon_s = (beta_fmincon_s(:, [1:end-1]) + beta_fmincon_s(:, [2:end]))/2;

fval_tau_qreg_s = X_mean_repeat'*beta_mid_qreg_s;
fval_tau_fmincon_s = X_mean_repeat'*beta_mid_fmincon_s;

[V_qreg_s,I_qreg_s] = sort(fval_tau_qreg_s);
[V_fmincon_s,I_fmincon_s] = sort(fval_tau_fmincon_s);

beta_sorted_qreg_s(:,[1:2:end]) = beta_qreg_s(:,I_qreg_s);
beta_sorted_qreg_s(:,[2:2:end]) = beta_qreg_s(:,I_qreg_s+1);
beta_sorted_fmincon_s(:,[1:2:end]) = beta_fmincon_s(:,I_fmincon_s);
beta_sorted_fmincon_s(:,[2:2:end]) = beta_fmincon_s(:,I_fmincon_s+1);


%% 4) Plot

for j_covar = 0:4
    
    figure; hold on;

    plot(taugrid_qreg,beta_qreg(j_covar+1,:),'g');
    plot(taugrid_break,beta_qreg_s(j_covar+1,:),'r');
   % plot(taugrid_break_sorted,beta_sorted_qreg_s(j_covar+1,:),'r--');
    
    plot(taugrid_break,beta_fmincon_s(j_covar+1,:),'b');
   % plot(taugrid_break_sorted,beta_sorted_fmincon_s(j_covar+1,:),'b--');
    
    % legend('Quantile Regression','QREG+fmincon', 'QREG+fmincon(Sorted)','fmincon+GA', 'fmincon+GA(Sorted)'); 
    legend('Quantile Regression','QREG+fmincon', 'fmincon+GA'); 
    title(sprintf('Data: %s. Beta %d',postfixdata, j_covar))
    
    print('-dpng','-r0',sprintf('beta%d_%s%s',j_covar,postfix,ntau_s));
end


