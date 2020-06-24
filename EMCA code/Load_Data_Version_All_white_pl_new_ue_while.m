% Piecewise linear MLE on white subsamples using 
%
% Updated 11/22/2018

clear;
close all;

% Part A) Set parameters
% data could be Angrist80 to Angrist00, or Jacob10
data=csvread('Angrist90.csv',1,0);
postfix = 'Angrist90';

% n_repeatsample could be 1
n_repeatsample = 1;
ntau = 15;
n_subsample = 1; % The number of observations in each subsample is calculated in Part B)
% Number of WLS iterations
n_WLS_iter = 200;

% taugrid is the old tau grid for quantile regression
% taugrid_midpoint is the mid point of tau segments
% taugrid_ue is the new uneven grid
% taugrid = (1:ntau)/(ntau+1);

taugrid_midpoint = (1/2*([0:(ntau-1)]  + [1:(ntau)]))/ntau;
taugrid_ue = taugrid_midpoint;
taugrid_ue(1) = 0;
taugrid_ue(end) = 1;


% Part B) Load original data & keep WHITE observations
y_orig_1 = data(:,1);
X_orig_1 = [ones(size(y_orig_1)), data(:,2:end)];
[nsample_orig_1, ncovar_orig_1] = size(X_orig_1);

% Keep WHITE observations
y_orig = y_orig_1(X_orig_1(:,ncovar_orig_1)==0);
X_orig = X_orig_1(X_orig_1(:,ncovar_orig_1)==0,:);
X_orig = X_orig(:,[1:(end-1)]);

% Calculate the number of observations in the subsamples
[nsample_orig, ncovar] = size(X_orig);
nsample = floor(1/n_subsample*nsample_orig);
if n_subsample == 17
    nsample = floor(500*nsample_orig^(2/5));
end
if n_subsample == 19
    nsample = floor(200*nsample_orig^(2/5));
end
% If n_subsample is 21, we do random sampling with replacement using the original number of observations
if n_subsample == 21
    nsample = nsample_orig;
end

% Part C) Preallocation of result variables
%number of mixture components, using mixture of normals
nmixtures=3;
%nvars is the total number of parameters
nvars=ncovar*ntau+3*nmixtures-2;

% qreg and WLS
recorder_qreg_repeat = nan(1,ncovar*ntau,n_repeatsample);
recorder_WLS_repeat = nan(1,ncovar*ntau,n_repeatsample);

% WLS result
betahatMatrix = nan(ncovar,ntau,n_WLS_iter,n_repeatsample);

% WLS start
recorder_WLS_start_repeat = nan(n_repeatsample,nvars);
fval_recorder_WLS_start_repeat = nan(1,n_repeatsample);
exit_recorder_WLS_start_repeat = nan(1,n_repeatsample);

% piecewise constant start to do piecewise linear
recorder_pc_start_repeat = nan(n_repeatsample,nvars);
fval_recorder_pc_start_repeat = nan(1,n_repeatsample);
exit_recorder_pc_start_repeat = nan(1,n_repeatsample);

% Mean of covariates
X_mean_repeat = nan(ncovar,n_repeatsample);

lower=[zeros(1,(ncovar*ntau))-10000, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(ncovar*ntau))+10000), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
lower_pl=[zeros(1,(ncovar*ntau))-10, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper_pl=[(zeros(1,(ncovar*ntau))+100), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];

for j = [1:n_repeatsample]   
    % Part D) Create a random subsample
    j
    options=optimoptions(@fmincon,'GradObj','on');
    opts = optimoptions('ga', 'MaxGenerations',500,'PopulationSize',500);
    
    if n_subsample ~= 1
        v_sample= randsample(nsample_orig,nsample,true);
    else
        v_sample = [1:nsample_orig];
    end
    [v_sample_sorted] = sort(v_sample);
    y = y_orig(v_sample);
    X = X_orig(v_sample,:);
    
    X_mean_repeat(:,j) = mean(X,1);
    
    % Part E) qreg and WLS
    Y = repmat(y, [1 ntau]);
    [fit] = quantlsfVector(X,y,taugrid_midpoint);
    fit_1 = reshape(fit, ncovar*ntau,1);
    recorder_qreg = fit_1';
    recorder_qreg_repeat(:,:,j) = recorder_qreg;
    
    % Preallocation for the matrix of probabilities
    pdfmatrix = nan(nsample,ntau);
    betahat   = fit';
    ehat = Y-X*betahat;
    sigmahat=std(y-X*quantlsf(X,y,.5)); % use only std of median regression residuals
    
    for j_WLS = 1 : 200
        pdfmatrix = normpdf(ehat./sigmahat);
        % the weight matrix is also nsample * ntau
        weight=ntau*pdfmatrix./repmat(sum(pdfmatrix,2),[1 ntau]);
        % j_WLS
        for j_tau = [1:ntau]
            weightV = weight(:,j_tau)';
            for j_covar = [1:ncovar]
                weightedX(:,j_covar) = weightV.*(X(:,j_covar)');
            end
            weightedy = weightV'.*y;
            betahat(:,j_tau)=(X'*weightedX)\(X'*weightedy);
        end
        betahatMatrix(:,:,j_WLS,j) = betahat;
        ehat=Y-X*betahat;
        sigmahat=(sum(sum(weight.*(ehat.^2)))/nsample/ntau)^.5;
        sigmahat_v(j_WLS,j) = sigmahat;
        if j_WLS > 1
        error_beta = norm(squeeze(betahatMatrix(:,:,j_WLS,j) - betahatMatrix(:,:,j_WLS-1,j)));
        if error_beta < 0.001
            break
        end
        end
    end
    recorder_WLS = reshape(betahat',ncovar*ntau,1)';
    recorder_WLS_repeat(:,:,j) = recorder_WLS;
    
    % Part F) MLE
  
    % Constants
    b=1;
    A = zeros(1,nvars);
    A(ncovar*ntau+1) = 1;
    A(ncovar*ntau+2) = 1;
    
    betahat = squeeze(betahatMatrix(:,:,j_WLS,j));
    recorder_WLS = reshape(betahat',ncovar*ntau,1)';
    start = [recorder_WLS, para_dist_default];
    % start = [fit_1',para_dist_default];
    [fit_hat,fval,exitflag] = fmincon(@(x)gradllfCovarparavector(x, ntau, nsample, nmixtures,1,y, X),start,A,b,[],[],lower,upper,[],options);
    
    recorder_WLS_start_repeat(j,:) = fit_hat;
    fval_recorder_WLS_start_repeat(j) = fval;
    exit_recorder_WLS_start_repeat(j) = exitflag;
    
    
    % G) Piecewise linear MLE
    % G-1) Sort piecewise constant result
    X_mean = mean(X);
    recorder_WLS_start_reshape = reshape(fit_hat(1:ntau*ncovar),ntau,ncovar);
    [beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,ncovar);
    beta_WLS_start_sorted = beta_WLS_start_sorted';
    
    beta_WLS_start_sorted(:, 1) = beta_WLS_start_sorted(:, 1) - 1/2*(beta_WLS_start_sorted(:, 2) - beta_WLS_start_sorted(:, 1));
    beta_WLS_start_sorted(:, end) = beta_WLS_start_sorted(:, end) + 1/2*(beta_WLS_start_sorted(:, end) - beta_WLS_start_sorted(:, end-1));
    
    % G-2) Construct the start value
    fit_1_temp = nan(1,ncovar*(ntau));
    for j_covar = [1:ncovar]
        % The start and the end index
        para_start = 2 + (j_covar-1)*(ntau);
        para_end = j_covar*(ntau);
        % The first constant
        fit_1_temp(1, para_start-1) = beta_WLS_start_sorted(j_covar, 1);
        % All increments after the first constant.
        fit_1_temp(1, [para_start : para_end]) = beta_WLS_start_sorted(j_covar, [2:end]) - beta_WLS_start_sorted(j_covar, [1:end-1]);
    end
    
    start = [fit_1_temp, fit_hat([(end - 3*nmixtures + 3) : end])];
    
    [fit_hat,fval,exitflag,output] = fmincon(@(x)gradl_CDF_Lei_GA_ue(x,  taugrid_ue, nmixtures, y', X'), start, A, b,[],[],lower_pl,upper_pl,[],options);
    recorder_pc_start_repeat(j,:) = fit_hat;
    fval_recorder_pc_start_repeat(j) = fval;
    exit_recorder_pc_start_repeat(j) = exitflag;

end
clear Y Wk weight pdfmatrix ehat A b



save (sprintf('RData_pl_%s_sub%d_ntau%d_white%d_ue_while',postfix,n_subsample,ntau,n_repeatsample))
