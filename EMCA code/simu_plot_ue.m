% This code plots the simulation results of piecewise linear MLE simulations
% Grids used: [1:ntau]/(ntau+1) for qreg, mid points for piecewise constant
% MLE, mid points extended to 0 and 1 for piecewise linear MLE.
% Haoyang Liu
% 3/17/2020
% Version note: because the Haas server was down, this will become the primary 
% plotting code 
% The following are informal notes while Haoyang edited the code on 3/17/2020
% 

clear
close all

% Dummy for if a case involves linear 
% 0 for PC only; 1 for PL only; 2 for PL plus WLS; 3 for y-star (a special simulation) 
pl = 0;

% EIV = "T PC"

% EIV = '3mix 99 knot PC'
% EIV = '3mix w WLS';
% EIV = '2mix PC'
% EIV = '3mix w WLS squared ystar';
 EIV = 'Laplace PC';
% EIV = 'Laplace PL';
% A) Load data
switch EIV
    case '3mix w WLS squared ystar'
        load Ye_GA_3mixS_WLS_ue_ystar
        pl = 3;
    case '3mix w WLS squared'
        load Ye_GA_3mixS_WLS_ue_1 % _1 is for uniformly distributed x1 
        %load Ye_GA_3mixS_WLS_ue % w/o _1 is for log normal x1 
    case '1 normal'
        load Ye_GA_1_WLS_ue 
        pl = 2;
        load Ye_GA_1_3_WLS_ue_qreg taugrid recorder_qreg_repeat_reshape
    case '3mix w WLS'
        load Ye_GA_3mix_3_WLS_ue_combined
        load Ye_GA_3mix_3_WLS_ue_99knot_500 taugrid recorder_qreg_repeat_reshape
        pl = 1;
        % case '3mix wo WLS'
        % Result from MLE_Ye_CDF_3mix_par_start_ue, which feeds qreg as the start value for piecewise constant MLE
        % load Ye_GA_3mix_3_ue
        % pl = 1;
    case '2mix PC'
        % Estimating a mixture 2 normals as a mixture of 3 normals.
        load Ye_GA_2mix_3_WLS_ue_500
        % Load an alternative qreg result because we need qreg to cover the
        % whole range of 0.01 to 0.99
        load Ye_GA_2mix_3_WLS_ue_qreg taugrid recorder_qreg_repeat_reshape
    case 'Laplace PC'
        % Estimating a mixture 2 normals as a mixture of 3 normals.
        load Ye_GA_L_3_WLS_ue_600
        recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted([1:500], :, :);
        load Ye_GA_L_3_WLS_ue_qreg_500 taugrid recorder_qreg_repeat_reshape
    case 'Laplace PL'
        % Estimating a Laplace as a mixture of 3 normals using piecewise
        % linear
        load Ye_GA_L_3_WLS_ue_pl_combined   
        pl = 1;
        load Ye_GA_L_3_WLS_ue_qreg taugrid recorder_qreg_repeat_reshape    
    case '3mix dummy PC'
        load Ye_GA_3mix_3_WLS_ue_dummy
        load Ye_GA_3mix_3_WLS_ue_dummy_qreg taugrid recorder_qreg_repeat_reshape
    case 'T PC'
        % Old data from the first round R&R
        % load Ye_GA_T_3_WLS_ue
        load Ye_GA_T_3_WLS_ue_combined
        recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted(exit_recorder_WLS_start_repeat ==  1, :, :);
        recorder_WLS_start_repeat_sorted = recorder_WLS_start_repeat_sorted([1:500], :, :);
        load Ye_GA_T_3_WLS_ue_qreg_500 taugrid recorder_qreg_repeat_reshape
    case '3mix 99 knot PC'      
        load Ye_GA_3mix_3_WLS_ue_99knot_500
    case '15 test CP'
        % I think this is the case for coverage probabilities for 15 knots
        load Ye_GA_3mix_3_par_start_3_cover_new_5
        [taugrid, taugrid_midpoint] = calculate_grid(ntau);
        taugrid_ue = taugrid_midpoint;
        recorder_qreg_repeat_reshape = nan(iter,ncovar,ntau);
        recorder_WLS_start_repeat_sorted = nan(iter,ncovar,ntau);
        
        for j_iter = [1:iter]
            recorder_qreg_repeat_reshape(j_iter,:,:) = reshape(squeeze(recorder_qreg(j_iter,:)), ntau, ncovar)';
            recorder_WLS_start_repeat_sorted(j_iter,:,:) = reshape(squeeze(recorder_qreg_start_sorted(j_iter,:)),ntau, ncovar)';
        end
        recorder_pc_start_repeat_recon = recorder_WLS_start_repeat_sorted;
        n_repeatsample = iter;
    otherwise
end

beta_true = [b0(taugrid_ue);b1(taugrid_ue);b2(taugrid_ue)];

% B) Means and stds of estiamtes
beta_qreg_full = squeeze(mean(recorder_qreg_repeat_reshape));
% Always collect piecewise constant means and stds
beta_WLS_start_sorted_full = squeeze(mean(recorder_WLS_start_repeat_sorted));
beta_WLS_start_sorted_std = squeeze(std(recorder_WLS_start_repeat_sorted));
% If pl == 1, collect piecewise linear 
if pl == 1
    beta_pc_s_full = squeeze(mean(recorder_pc_start_repeat_recon));
    beta_pc_s_std = squeeze(std(recorder_pc_start_repeat_recon));
end
% if pl == 2, also collect WLS
if pl == 2
    beta_WLS_sorted_full = squeeze(mean(recorder_WLS_repeat_sorted));
	beta_WLS_sorted_std = squeeze(std(recorder_WLS_repeat_sorted));
end
% For the y-star case
if pl == 3
    beta_qreg_full_1 = squeeze(mean(recorder_qreg_repeat_reshape_1));
end

% C) Plot
for j_covar = [0:ncovar-1]
    figure;
    hold on;
    
    % Always need to plot the truth and the qreg estimate 
    plot(taugrid_ue,beta_true(j_covar+1,:),'b');   
    plot(taugrid,beta_qreg_full(j_covar+1,:),'m','LineWidth',1);
    
    % For the special y-star case, we also need to plot quantile regression
    % for y-star
    if pl == 3
        plot(taugrid,beta_qreg_full_1(j_covar+1,:),'m--','LineWidth',1); 
    end
    
    % For the piecewise linear case, plot piecewise linear. 
    % Otherwise, always plot piecewise constant 
    if pl == 1
        plot(taugrid_ue,beta_pc_s_full(j_covar+1,:),'r','LineWidth',1);
    else
        plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:),'k','LineWidth',1);
    end
	
    % If pl == 2, plot WLS
	if pl == 2
		plot(taugrid_midpoint,beta_WLS_sorted_full(j_covar+1,:),'m','LineWidth',1);
	end
    
    % Plot confidence intervals
    if pl == 1
        plot(taugrid_ue,beta_pc_s_full(j_covar+1,:) + 1.96*beta_pc_s_std(j_covar+1,:),'r--');
        plot(taugrid_ue,beta_pc_s_full(j_covar+1,:) - 1.96*beta_pc_s_std(j_covar+1,:),'r--');
    else
        plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:) + 1.96*beta_WLS_start_sorted_std(j_covar+1,:),'k--');
        plot(taugrid_midpoint,beta_WLS_start_sorted_full(j_covar+1,:) - 1.96*beta_WLS_start_sorted_std(j_covar+1,:),'k--');
    end	
	if pl == 2
       plot(taugrid_midpoint,beta_WLS_sorted_full(j_covar+1,:) + 1.96*beta_WLS_sorted_std(j_covar+1,:),'m--');
       plot(taugrid_midpoint,beta_WLS_sorted_full(j_covar+1,:) - 1.96*beta_WLS_sorted_std(j_covar+1,:),'m--');	
	end
	
	switch pl 
        case 0
		legend('truth','quantile regression','piecewise constant');
    	case 1
		legend('truth','quantile regression','piecewise linear(fmincon)');
		case 2
		legend('truth','quantile regression','piecewise constant(fmincon)','WLS');
	    case 3 
        legend('truth','quantile regression','quantile regression ystar','piecewise linear(fmincon)');    
	end
    
    title(sprintf('EIV: %s. # of simulations: %d',EIV, min(n_repeatsample,j)))
    print('-dpng','-r0',sprintf('Sim_%s_%d_%d',EIV,min(n_repeatsample,j),j_covar));
end

beta_WLS_start_sorted_full(:,1) = beta_WLS_start_sorted_full(:,1) - 1/2*(beta_WLS_start_sorted_full(:,2) - beta_WLS_start_sorted_full(:,1));
beta_WLS_start_sorted_full(:,end) = beta_WLS_start_sorted_full(:, end) + 1/2*(beta_WLS_start_sorted_full(:, end) - beta_WLS_start_sorted_full(:, end-1));

beta_WLS_start_sorted_std(:,1) = beta_WLS_start_sorted_std(:,1) - 1/2*(beta_WLS_start_sorted_std(:,2) - beta_WLS_start_sorted_std(:,1));
beta_WLS_start_sorted_std(:,end) = beta_WLS_start_sorted_std(:, end) + 1/2*(beta_WLS_start_sorted_std(:, end) - beta_WLS_start_sorted_std(:, end-1));


if pl == 1
    output_matrix1 = [taugrid_ue;beta_true; beta_pc_s_full; beta_pc_s_full(:,:) + 1.96*beta_pc_s_std(:,:); beta_pc_s_full(:,:) - 1.96*beta_pc_s_std(:,:)]';
else
    output_matrix1 = [taugrid_ue;beta_true; beta_WLS_start_sorted_full; beta_WLS_start_sorted_full(:,:) + 1.96*beta_WLS_start_sorted_std(:,:); beta_WLS_start_sorted_full(:,:) - 1.96*beta_WLS_start_sorted_std(:,:)]';
end


output_matrix2 = [taugrid; beta_qreg_full]';

csvwrite(sprintf('result_%s_1.csv',EIV), output_matrix1);
csvwrite(sprintf('result_%s_2.csv',EIV), output_matrix2);


if pl == 2
    output_matrix3 = [taugrid_midpoint;beta_WLS_sorted_full; beta_WLS_sorted_full(:,:) + 1.96*beta_WLS_sorted_std(:,:); beta_WLS_sorted_full(:,:) - 1.96*beta_WLS_sorted_std(:,:)]';
    csvwrite(sprintf('result_%s_3.csv',EIV), output_matrix3);
end


% This part calculate 


if pl == 3
    x1vector = [0:0.2:2*2.34];
    x2vector = 2*2.34*(x1vector/(2*2.34)).^2;
  
    for j_tau = [2 8 14]
        figure; 
        hold on; 
        plot(x1vector, beta_true(1,j_tau) + beta_true(2,j_tau)*x1vector + beta_true(3,j_tau)*x2vector,'b');
        plot(x1vector, beta_qreg_full(1,j_tau) + beta_qreg_full(2,j_tau)*x1vector,'m');
        plot(x1vector, beta_qreg_full_1(1,j_tau) + beta_qreg_full_1(2,j_tau)*x1vector,'m--'); 
        plot(x1vector, beta_WLS_start_sorted_full(1,j_tau) + beta_WLS_start_sorted_full(2,j_tau)*x1vector,'k'); 
        
       lgd = legend('truth','quantile regression','quantile regression ystar','piecewise linear(fmincon)');  
       c = lgd.Location;
       lgd.Location = 'northwest';
       title(sprintf('Q-tau(Y|X) tau = %0.1g',taugrid_ue(j_tau)))
       print('-dpng','-r200',sprintf('qyx%d',j_tau));
 
    end
 
end
