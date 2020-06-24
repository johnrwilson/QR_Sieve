% This codes implements Ye's approach on a mixture of 3 EIV
% It also specifies the qreg results as the start value for GA
% Haoyang Liu
% 4/21/2018

clear all;
close all;


% % %% 1) Definition of constants

opts = optimoptions('ga', 'MaxGenerations',5000,'PopulationSize',500,'UseParallel',true,'display','iter');
      
% Number of grid intervals. The number of grid points is one more than ntau
ntau = 10;

% In MLE, taugrid always covers the entire [0,1] interval. Thus the last
% grid point (1) is always omitted
taugrid = (0:(ntau))/(ntau);
taugrid_qreg = (0:(ntau))/(ntau);
epsilon = 0.01;
taugrid_qreg(1) = epsilon + taugrid_qreg(1);
taugrid_qreg(end) =  -epsilon + taugrid_qreg(end);

% Number of iterations (or runs of simulations)
iter = 10;
nsample = 100000;

%number of mixture components, using mixture of normals
nmixtures=3;
%number of true components
nmixtures_truth = 3;

% nvars is the total number of parameters
% 3 for cs (constants)
% 3*ntau for beta (there are ntau intervals, and thus ntau slopes to estimate
% 3*nmixtures-2 for the distributional parameters
nvars=3+3*ntau+3*nmixtures-2;

% Preallocation for both quantile regression and MLE
% For qreg, tau grid is from 0 to 1. Thus there are ntau+1 points in the grid
recorder_qreg = nan(iter,3*ntau+3);
% For qreg start, the number of parameters is calculated above
recorder_qreg_start = nan(iter,nvars);

x1matrix = nan(iter,nsample);
x2matrix = nan(iter,nsample);
ytmatrix = nan(iter,nsample);


% Constants
b=1;
A = zeros(1,nvars);
A(3*ntau+3+1) = 1;
A(3*ntau+3+2) = 1;


% Mixing probability of each component
lambda1_true = 0.5;
lambda2_true = 0.25;

% Mean of each component
mu1_true = -3; 
mu2_true = 2;

% Preprocess the parameters
[lambdavector_true,muvector_true,q_true] =  preprocesslambdamu([lambda1_true,lambda2_true],[mu1_true,mu2_true]);
lambda3_true = lambdavector_true(nmixtures_truth);
mu3_true = muvector_true(nmixtures_truth);

%true parameters for the error term:
sv_true=[lambda1_true,lambda2_true,mu1_true,mu2_true,1,1,1];


%Define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower = [-0.1,zeros(1,ntau),0.5,repmat(0.01,1,ntau),-0.1,repmat(0.01,1,ntau),(zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper = [1,repmat(0.1,1,ntau),2,ones(1,ntau),2,repmat(0.5,1,ntau),(zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];

for j_iter = [1:iter]
    j_iter
%% 2) Simulate the data
    % Betas are simulated in a continuous range. So no need to change that. 
	tau_simu = rand(1,nsample);         
    beta0_simu = b0(tau_simu);
	beta1_simu = b1(tau_simu);
	beta2_simu = b2(tau_simu);
    
	x1r = exp(randn(1,nsample)); 
    x2r = exp(randn(1,nsample));
    
    x1matrix(j_iter,:) = x1r;
    x2matrix(j_iter,:) = x2r;
    
    y_ntemp = zeros(1,nsample);
    j_sample_begin = 1;
    for j_mixture = [1 : nmixtures_truth]
        if j_mixture == nmixtures_truth
            j_sample_end = nsample;
        else
            j_sample_end = min( ceil(sum(lambdavector_true(1:j_mixture))*nsample) , nsample);
        end
        y_ntemp([ j_sample_begin:j_sample_end]) = normrnd(muvector_true(j_mixture),1,1,[(j_sample_end - j_sample_begin +1)]);
        j_sample_begin = j_sample_end+1;
    end
    [temp , y_n_index] = sort(rand(1,nsample));


    y_n = y_ntemp(y_n_index);
    
    
    y_s = beta0_simu + beta1_simu.*x1r + beta2_simu.*x2r;
    y = y_n+y_s;
    ytmatrix(j_iter,:) = y;        
    y = y';
    
%% 3) QREG
    X = [ones(1,nsample); x1r; x2r]';
    Y = repmat(y, [1 ntau]);
    
    [fit] = quantlsfVector(X,y,taugrid_qreg);        
    fit_1 = [squeeze(fit(:,1));squeeze(fit(:,2));squeeze(fit(:,3))];
    recorder_qreg(j_iter,:)=fit_1';
    
    
%% 4) MLE
      %
      
      fit_1_temp = fit_1';
      ncovar = 3;
      ntau_plus = ntau + 1;
      for j_covar = [1:ncovar]
          para_start = 2 + (j_covar-1)*ntau_plus;
          para_end = j_covar*ntau_plus;
          fit_1_temp(1, [para_start : para_end]) = fit_1_temp(1, [para_start : para_end]) - fit_1_temp(1, [para_start-1 : para_end-1]);       
      end
    
    
      para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];
      start = [fit_1_temp, para_dist_default];   
      
      opts.InitialPopulationMatrix = start;
      
      
      display('reached MLE')  
      [fit_hat,fval,exitflag,output] = ga(@(x)gradl_CDF_GA(x,  taugrid, nmixtures, y', X'), nvars, A, b, [],[],lower,upper,[],[],opts);
      recorder_qreg_start(j_iter,:) = fit_hat;
      fval_recorder_qreg_start(j_iter) = fval;
      exit_recorder_qreg_start(j_iter) = exitflag;
      iteration_recorder_qreg_start(j_iter) = output.generations;
      funcCount_recorder_qreg_start(j_iter) = output.funccount;
      
      
      
    
end
clear X Y

%% 5) Save the result
save Ye_GA_3mix_3_par_start
return;
%iter = 1
%% 6) Plot the result
a1rec = zeros(iter, ntau+1);
a2rec = zeros(iter, ntau+1);
a3rec = zeros(iter, ntau+1);

for j = 1:iter
    [b1,b2,b3,sigma] = reconstruct_beta2(recorder_qreg_start(j,[1:(3*ntau+3+1)]));
    a1rec(j,:) = b1;
    a2rec(j,:) = b2;
    a3rec(j,:) = b3;
end

ma1 = mean(a1rec,1);
ma2 = mean(a2rec,1);
ma3 = mean(a3rec,1);

figure;
plot(taugrid,ma1,'k-o',taugrid,zeros(1,length(taugrid)),'r-o',taugrid_qreg,mean(recorder_qreg([1:iter],1:11),1),'b-o');

legend('ours','true','q-reg','Location','northwest');
print('m1_3mix_3','-dpng');
figure;
plot(taugrid,ma2,'k-o',taugrid,exp(taugrid),'r-o',taugrid_qreg,mean(recorder_qreg([1:iter],12:22),1),'b-o');
%ylim([0.8 3]);
legend('ours','true','q-reg','Location','northwest');
print('m2_3mix_3','-dpng');
figure;
plot(taugrid,ma3,'k-o',taugrid,sqrt(taugrid),'r-o',taugrid_qreg,mean(recorder_qreg([1:iter],23:33),1),'b-o');
%ylim([-0.2 1.2]);
legend('ours','true','q-reg','Location','northwest');
print('m3_3mix_3','-dpng');
