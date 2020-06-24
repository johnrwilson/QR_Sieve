% This is the plotting code for GA real data results

% It plots sorted estimates
% It doesn't let sorting done at a finer grid 
% 3/25/2018
% Haoyang Liu

clear
close all 

%%% Part A) Set Parameters and Load Data
postfix = 'Angrist80';
error = 'Mixture of 3';
ntau = 33;
load (sprintf('Lei_GA_%s_ntau%d',postfix,ntau))
postfixdata = 'Angrist 80'

%%%% The following lines are for white 
% load GA_Angrist90_sub10000_ntau33_white
% postfix = 'Angrist90white';
% postfixdata = 'Angrist 90 White'

%load GA_Angrist90_sub3000_ntau33
%postfix = 'Angrist9030000';
%postfixdata = 'Angrist 90 3000';
%%% End of Part A)


%%% Part B) Reconstruct beta
beta = reconstruct_beta( (reshape(recorder_qreg_start_repeat(1,[1:(ncovar*ntau+ncovar)]), [ntau+1, ncovar]))'  );
beta_sorted = nan(ncovar, 2*ntau);

%%% Part C) Create "taugrid_break_sorted", which repeats each middle entry twice 
step_original = taugrid_qreg(3) - taugrid_qreg(2);
% taugrid_break is still a little different from taugrid_qreg because the
% first entry is 0 instead of 0+0.00001
taugrid_break = [0:step_original:1];
taugrid_break_sorted = repeat_HL(taugrid_break,2);
% Repeating it twice gives us end points
taugrid_break_sorted = taugrid_break_sorted(1:end-1);


%%% Part D) Fill in the parametes for qreg
beta_qreg = (reshape(recorder_qreg(1,[1:(ncovar*ntau+ncovar)]), [ntau+1, ncovar]))';

%%% Part E) Sorting 
beta_mid = (beta(:, [1:end-1]) + beta(:, [2:end]))/2;
fval_tau_qreg_start = X_mean_repeat'*beta_mid;
    

[V,I_qreg_start] = sort(fval_tau_qreg_start);
    
beta_sorted(:,[1:2:end]) = beta(:,I_qreg_start);
beta_sorted(:,[2:2:end]) = beta(:,I_qreg_start+1);

%% 4) Plot
figure; hold on;


plot(taugrid_qreg,beta_qreg(2,:),'g');
plot(taugrid_break,beta(2,:),'r');
plot(taugrid_break_sorted,beta_sorted(2,:),'r--');

legend('Quantile Regression','Piecewise Linear MLE', 'Piecewise Linear MLE (Sorted)'); 
title(sprintf('Data: %s',postfixdata))
%return;
print('-dpng','-r0',sprintf('beta1_Lei_GA_%s',postfix));


