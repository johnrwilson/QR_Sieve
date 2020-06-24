% This code plots the simulation results of piecewise linear MLE
% Simulations
% Haoyang Liu
% 8/14/2018

clear
close all

% A) Load data
load Ye_GA_T_3_par_start_combined_1 ; % Ye_GA_2mix_3_par_start Ye_GA_L_3_par_start_combined Ye_GA_3mix_3_par_start_combined Ye_GA_3mix_3_par_start_3_dummy
EIV = "T EIV"; %  "Two-normal EIV"  "Laplacian EIV" "Three-normal EIV" "Three-normal EIV (x1 is a dummy)"
eiv_output =  "T"; % "2mix" "L" "3mix" "3dummy"

% B) Number of iterations
iter = length(exit_recorder_pc_start);

% C) Create the level of the coefficients from the piecewise linear
% coefficients; Sort the piecewise constant coefficients
a1rec = zeros(iter, ntau);
a2rec = zeros(iter, ntau);
a3rec = zeros(iter, ntau);

for j = 1:iter
    [b1,b2,b3] = reconstruct_beta2(recorder_pc_start(j,[1:(3*ntau)]));
    a1rec(j,:) = b1;
    a2rec(j,:) = b2;
    a3rec(j,:) = b3;
  
    x1r = squeeze(x1matrix(j,:));
    x2r = squeeze(x2matrix(j,:));  
    X_mean = [1,mean(x1r),mean(x2r)];
    
    beta_qreg_start =  reshape(recorder_qreg_start(j,[1 : 3*ntau])',ntau,3);
    fval_tau_qreg_start = X_mean*beta_qreg_start';
    [V,I_qreg_start] = sort(fval_tau_qreg_start);
    for j_x = [1:3]
        beta_qreg_start_sorted(j,[(j_x-1)*ntau+1 : j_x*ntau]) = beta_qreg_start(I_qreg_start,j_x)';
    end
end
%a1rectemp = a1rec;
%a1rec = a1rec(a1rectemp(:,end)<5,:);
%a2rec = a2rec(a1rectemp(:,end)<5,:);
%a3rec = a3rec(a1rectemp(:,end)<5,:);

% F) Calculate averages and standard deviations
ma1_pl = mean(a1rec,1);
ma2_pl = mean(a2rec,1);
ma3_pl = mean(a3rec,1);
ma1qreg = mean(recorder_qreg([1:iter],[1:(ntau)]),1);
ma2qreg = mean(recorder_qreg([1:iter],[ntau+1:(2*ntau)]),1);
ma3qreg = mean(recorder_qreg([1:iter],[2*ntau+1:(3*ntau)]),1);
ma1_pc = mean(beta_qreg_start_sorted([1:iter],[1:(ntau)]),1);
ma2_pc = mean(beta_qreg_start_sorted([1:iter],[ntau+1:(2*ntau)]),1);
ma3_pc = mean(beta_qreg_start_sorted([1:iter],[2*ntau+1:(3*ntau)]),1);


std1_pl = std(a1rec);
std2_pl = std(a2rec);
std3_pl = std(a3rec);
std1qreg = std(recorder_qreg([1:iter],[1:(ntau)]),1);
std2qreg = std(recorder_qreg([1:iter],[ntau+1:(2*ntau)]),1);
std3qreg = std(recorder_qreg([1:iter],[2*ntau+1:(3*ntau)]),1);
std1_pc = std(beta_qreg_start_sorted([1:iter],[1:(ntau)]),1);
std2_pc = std(beta_qreg_start_sorted([1:iter],[ntau+1:(2*ntau)]),1);
std3_pc = std(beta_qreg_start_sorted([1:iter],[2*ntau+1:(3*ntau)]),1);

% Preallocate memory for the MSE variable. 
% The first 3 is for the number of estimators. The second 3 is for the three betas. The third dimension is for the number of taus 
MSE = zeros(ntau,7);
MSE(:,1) = taugrid';
MSE(:,2) = (std1qreg.^2 + (ma1qreg - beta0_true).^2)';
MSE(:,3) = (std1_pl.^2 + (ma1_pl - beta0_true).^2)';
MSE(:,4) = (std2qreg.^2 + (ma2qreg - beta1_true).^2)';
MSE(:,5) = (std2_pl.^2 + (ma2_pl - beta1_true).^2)';
MSE(:,6) = (std3qreg.^2 + (ma3qreg - beta2_true).^2)';
MSE(:,7) = (std3_pl.^2 + (ma3_pl - beta2_true).^2)';

csvwrite(sprintf('MSE%s.csv',eiv_output),MSE,2,0);


% The following lines are temporarily commented because we already
% produced those graphs. We don't need to redo the graphs for now
figure;
hold on;
plot(taugrid,ma1_pl,'k',taugrid,ma1_pc,'b',taugrid,ma1qreg,'m',taugrid,beta0_true,'r');
plot(taugrid, ma1_pl-1.96*std1_pl,'k--')
plot(taugrid, ma1_pl+1.96*std1_pl,'k--')
plot(taugrid, ma1_pc-1.96*std1_pc,'b--')
plot(taugrid, ma1_pc+1.96*std1_pc,'b--')
legend('PL MLE','PC MLE','q-reg','truth','Location','northwest');
axis([0.05 0.95 min(ma1qreg)-0.5 max(ma1qreg)+0.5])
title(sprintf('beta 0, %d simulations %s',iter, EIV))
print(sprintf('m0_%s_%d',eiv_output,iter),'-dpng');
results(1,:) = taugrid;
results(2,:) = beta0_true;
results(3,:) = ma1_pl;
results(4,:) = ma1_pc;
results(5,:) = ma1qreg;
results(6,:) = ma1_pl-1.96*std1_pl;
results(7,:) = ma1_pl+1.96*std1_pl;
results(8,:) = ma1_pc-1.96*std1_pc;
results(9,:) = ma1_pc+1.96*std1_pc;




figure;
hold on;
plot(taugrid,ma2_pl,'k',taugrid,ma2_pc,'b',taugrid,ma2qreg,'m',taugrid,beta1_true,'r');
plot(taugrid, ma2_pl-1.96*std2_pl,'k--')
plot(taugrid, ma2_pl+1.96*std2_pl,'k--')
plot(taugrid, ma2_pc-1.96*std2_pc,'b--')
plot(taugrid, ma2_pc+1.96*std2_pc,'b--')

legend('PL MLE','PC MLE','q-reg','truth','Location','northwest');
axis([0.05 0.95 min(beta1_true)-0.5 max(beta1_true)+0.5])
title(sprintf('beta 1, %d simulations  %s',iter, EIV))
print(sprintf('m1_%s_%d',eiv_output,iter),'-dpng');

results(10,:) = beta1_true;
results(11,:) = ma2_pl;
results(12,:) = ma2_pc;
results(13,:) = ma2qreg;
results(14,:) = ma2_pl-1.96*std2_pl;
results(15,:) = ma2_pl+1.96*std2_pl;
results(16,:) = ma2_pc-1.96*std2_pc;
results(17,:) = ma2_pc+1.96*std2_pc;


figure;
hold on;
plot(taugrid,ma3_pl,'k',taugrid,ma3_pc,'b',taugrid,ma3qreg,'m',taugrid,beta2_true,'r');
plot(taugrid, ma3_pl-1.96*std3_pl,'k--')
plot(taugrid, ma3_pl+1.96*std3_pl,'k--')
plot(taugrid, ma3_pc-1.96*std3_pc,'b--')
plot(taugrid, ma3_pc+1.96*std3_pc,'b--')
%plot(taugrid, ma3qreg-1.96*std3qreg,'m--')
%plot(taugrid, ma3qreg+1.96*std3qreg,'m--')
legend('PL MLE','PC MLE','q-reg','truth','Location','northwest');
axis([0.05 0.95 min(beta2_true)-0.5 max(beta2_true)+0.5])
title(sprintf('beta 2, %d simulations  %s',iter, EIV))
print(sprintf('m2_%s_%d',eiv_output,iter),'-dpng');


results(18,:) = beta2_true;
results(19,:) = ma3_pl;
results(20,:) = ma3_pc;
results(21,:) = ma3qreg;
results(22,:) = ma3_pl-1.96*std3_pl;
results(23,:) = ma3_pl+1.96*std3_pl;
results(24,:) = ma3_pc-1.96*std3_pc;
results(25,:) = ma3_pc+1.96*std3_pc;

results = results';

% The following line is commented to not overwrite existing files (already
% with the right headers)
csvwrite(sprintf('%s.csv',eiv_output),results,2,0);



