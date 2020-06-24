% This code plots the piecewise LINEAR full sample real data result on an
% uneven grid for MLE, an even grid for quantile regression


% Updated on 12/9/2018
% Haoyang Liu

clear
close all 

  %%% Part A) Set Parameters and load Data
    postfix =  'Jacob10';  %  'Angrist00' %  
    postfixdata =   'Jacob 10'; %  'Angrist 00' %  

    switch postfix
        case 'Angrist80'
            load RData_pl_Angrist80_sub1_ntau15_white1_ue1
        case 'Angrist90'
            load RData_pl_Angrist90_sub1_ntau15_white1_ue1
        case 'Angrist00'
            load RData_pl_Angrist00_sub1_ntau15_white1_ue1
        otherwise
            load RData_pl_Jacob10_sub1_ntau15_white1_ue1
    end




  %%% Part B) Plot the estimates
    for j_covar = [1]    
        figure; hold on;
    
        plot(taugrid,squeeze(fit(:,j_covar+1)),'g');
        plot(taugrid_midpoint,beta_WLS_start_sorted(j_covar+1,:),'k');
        plot(taugrid_ue,beta_pc_s(j_covar+1,:),'r');
     
    
        legend('quantile regression','piecewise constant','piecewise linear(fmincon)','Location','northeast');% , 'piecewise linear(GA)'
        title(sprintf('Data: %s. Number of knots %d',postfixdata, (ntau)))
        print('-dpng','-r0',sprintf('beta%d_%s_ue',(ntau),postfix));
    end
    
   % save pl_Jacob10_sub1_ntau15_white1_result taugrid_midpoint taugrid_ue  beta_qreg beta_WLS_start_sorted beta_pc_s postfixdata postfix


