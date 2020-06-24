% Matlab translation of R code deoptimpar_cjp.r
% This version first uses regular quantile regression and then do WLS. 
% The result of WLS is used as the starting point for MLE
% cd "bulk/qreg-eiv/Codes/Current version"
% 
% Haoyang Liu
% 1/15/2014

clear all;
close all;

% Call for multiple workers



do_MLE_WLS_start = 0;
do_MLE_qreg_start = 1;
do_MLE_truth_start = 0;

%% 1) Definition of constants
% Number of grid points
ntau = 9;
% taugrid is a grid on [0,1]
taugrid = (1:ntau)/(ntau+1);

% Number of iterations (or runs of simulations)
iter = 200;

nsample = 10000;

% Number of WLS iterations
n_WLS_iter = 80;

%number of mixture components, using mixture of normals
nmixtures=3;
nmixtures_truth = 3;
%nvars is the total number of parameters
nvars=3*ntau+3*nmixtures-2;

% Option for fmincon
options=optimoptions(@fmincon,'display','iter','GradObj','on','UseParallel',true);%,'DerivativeCheck','on'

% Default parameter for the distributions of measurement errors. They are
% used together with the result of WLS as the starting point for MLE
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];
% Preallocation
recorder_qreg = nan(iter,3*ntau);
recorder_WLS = nan(iter,3*ntau);

recorder_ParaDist = nan(iter,(nvars-3*ntau));
if do_MLE_WLS_start == 1
    recorder_WLS_start = nan(iter,nvars);
    fval_recorder_WLS_start = nan(iter,1);
    exit_recorder_WLS_start = nan(iter,1);
    iteration_recorder_WLS_start = nan(iter,1);
    funcCount_recorder_WLS_start = nan(iter,1);
    firstorderopt_recorder_WLS_start = nan(iter,1);
end

if do_MLE_qreg_start == 1
    recorder_qreg_start = nan(iter,nvars);
    fval_recorder_qreg_start = nan(iter,1);
    exit_recorder_qreg_start = nan(iter,1);
    iteration_recorder_qreg_start = nan(iter,1);
    funcCount_recorder_qreg_start = nan(iter,1);
    firstorderopt_recorder_qreg_start = nan(iter,1);
end

if do_MLE_truth_start == 1
    recorder_truth_start = nan(iter,nvars);
    fval_recorder_truth_start = nan(iter,1);
    exit_recorder_truth_start = nan(iter,1);
    iteration_recorder_truth_start = nan(iter,1);
    funcCount_recorder_truth_start = nan(iter,1);
    firstorderopt_recorder_truth_start = nan(iter,1);
end



x1matrix = nan(iter,nsample);
x2matrix = nan(iter,nsample);
ytmatrix = nan(iter,nsample);

% Constants
b=1;
A = zeros(1,nvars);
A(3*ntau+1) = 1;
A(3*ntau+2) = 1;
A(3*ntau+3) = 0;

b_ParaDist = 1;
A_ParaDist = zeros(1,(nvars-3*ntau));
A_ParaDist(1) = 1;
A_ParaDist(2) = 1;
A_ParaDist(3) = 0;


% true coefficients
beta0_true = b0(taugrid);
beta1_true = b1(taugrid);
beta2_true = b2(taugrid);

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

%true parameter sv:
sv_true=[lambda1_true,lambda2_true,mu1_true,mu2_true,1,1,1];

%define some rough lower and upper bounds for (beta,sigma). minimum weight of component=0.01, maximum=1, minimum mean=-10,maximum=10,
%minimum st.d=0.01,maximum=10
lower=[zeros(1,(3*ntau)), (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(3*ntau))+3), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
lower_ParaDist = [(zeros(1,(nmixtures-1))+0.01),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
upper_ParaDist = [(zeros(1,(nmixtures-1))*0+0.97),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];

for j_iter = [1:iter]
    display('got here')
%% 2) Simulate the data
	tau_simu = rand(1,nsample);         
    beta0_simu = b0(tau_simu);
	beta1_simu = b1(tau_simu);
	beta2_simu = b2(tau_simu);
    
	x1r = exp(normrnd(0,1,1,nsample));
    x2r = exp(normrnd(0,1,1,nsample));
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
    
    %% 3) WLS
    X = [ones(1,nsample); x1r; x2r]';
    Y = repmat(y, [1 ntau]);
    
    [fit] = quantlsfVector(X,y,taugrid);        
    fit_1 = [squeeze(fit(:,1));squeeze(fit(:,2));squeeze(fit(:,3))];
    recorder_qreg(j_iter,:)=fit_1';
    
    % Preallocation for the matrix of probabilities
    pdfmatrix = nan(nsample,ntau);    
    
    betahat   = fit';
    ehat = Y-X*betahat;    
    
  % sigmahat=std(y-X*quantlsf(X,y,.5)); % use only std of median regression residuals
   sigmahat=(sum(sum((ehat.^2)))/nsample/ntau)^.5; % an alternative to try.
    
    for j_WLS=1:n_WLS_iter
        % the weight matrix is also nsample * ntau
        weight=ntau*normpdf(ehat./sigmahat)./repmat(sum(normpdf(ehat./sigmahat),2),[1 ntau]);
        % j_WLS
        for j_tau = [1:ntau]
            Wk=sparse(diag(weight(:,j_tau)));
            betahat(:,j_tau)=(X'*Wk*X)\(X'*Wk*y);          
        end
        ehat=Y-X*betahat;
		sigmahat=(sum(sum(weight.*(ehat.^2)))/nsample/ntau)^.5;
    end
    recorder_WLS(j_iter,:) = [squeeze(betahat(1,:)),squeeze(betahat(2,:)),squeeze(betahat(3,:))];
    
    ehat1 = y-X*quantlsf(X,y,.5);
   % [ParaDist_hat] = fmincon(@(x)gradllf_WLS_resid(x, nsample, nmixtures, 1, ehat1),para_dist_default,A_ParaDist,b_ParaDist,[],[],lower_ParaDist,upper_ParaDist,[],options);
    
   % recorder_ParaDist(j_iter,:) = ParaDist_hat;
     
%     %% 4) MLE using WLS as start
%     if do_MLE_WLS_start ==1
%         start = [squeeze(recorder_WLS(j_iter,:)),para_dist_default];
%          [ start ] = changestart( start,lower,upper );
%         if gradllf(start, ntau, nsample, nmixtures, y, x1r, x2r) == Inf
%             exit_recorder_WLS_start(j_iter) = -10;
%         else
%         [fit_hat,fval,exitflag,output] = fmincon(@(x)gradllf(x, ntau, nsample, nmixtures, y, x1r, x2r),start,A,b,[],[],lower,upper,[],options);
%         j_iter
%         
%         recorder_WLS_start(j_iter,:)=fit_hat;
%         fval_recorder_WLS_start(j_iter) = fval;
%         exit_recorder_WLS_start(j_iter) = exitflag;
%         iteration_recorder_WLS_start(j_iter) = output.iterations;
%         funcCount_recorder_WLS_start(j_iter)= output.funcCount;
%         firstorderopt_recorder_WLS_start(j_iter) = output.firstorderopt;
%         end
%     end

    
 
    %% 5) MLE using qreg as start
    if do_MLE_qreg_start == 1
        start = [squeeze(recorder_qreg(j_iter,:)),para_dist_default];
         [ start ] = changestart( start,lower,upper );
        if gradllfCovarparavector(start, ntau, nsample, nmixtures, y, x1r, x2r) == Inf
            exit_recorder_qreg_start(j_iter) = -10;
        else
                xmatrix = [ones(1,nsample);x1r;x2r]';
        [fit_hat,fval,exitflag,output] = fmincon(@(x)gradllfCovarparavector(x,  ntau, nsample, nmixtures,1, y, xmatrix),start,A,b,[],[],lower,upper,[],options);
     
        
        recorder_qreg_start(j_iter,:)=fit_hat;
        fval_recorder_qreg_start(j_iter) = fval;
        exit_recorder_qreg_start(j_iter) = exitflag;
        iteration_recorder_qreg_start(j_iter) = output.iterations;
        funcCount_recorder_qreg_start(j_iter) = output.funcCount;
        firstorderopt_recorder_qreg_start(j_iter) = output.firstorderopt;
        end
    end
    
    %% 6) MLE using the truth as start
%     if do_MLE_truth_start == 1
%         start = Para_true;
%         
%         [fit_hat,fval,exitflag,output] = fmincon(@(x)gradllf(x, ntau, nsample, nmixtures, y, x1r, x2r),start,A,b,[],[],lower,upper,[],options);
%         
%         
%         recorder_truth_start(j_iter,:)=fit_hat;
%         fval_recorder_truth_start(j_iter) = fval;
%         exit_recorder_truth_start(j_iter) = exitflag;
%         iteration_recorder_truth_start(j_iter) = output.iterations;
%         funcCount_recorder_truth_start(j_iter) = output.funcCount;
%         firstorderopt_recorder_truth_start(j_iter) = output.firstorderopt;
%     end
    
end

save Test_MLE_Residuals_Median_three_as3_logn_E_Vector

