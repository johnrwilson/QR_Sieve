% Piecewise linear MLE on white subsamples using 
% This is the final version using a grid like [0.1:0.1:0.9] for quantile
% regression, and connecting the mid point values for pieceswise constant to 
% construct the start values for piecewise linear MLE
% Updated 12/9/2018

clear;
close all;

% Part A) Set parameters
% data could be Angrist80 to Angrist00, or Jacob10
data=csvread('Angrist00.csv',1,0);
postfix = 'Angrist00';

% n_repeatsample could be 1
n_repeatsample = 500;
ntau = 15;
n_subsample = 21; % The number of observations in each subsample is calculated in Part B)
% When n_subsample == 21, it runs the bootstrap procedure with the number
% of observations same as in the full sample. 

% Number of WLS iterations
n_WLS_iter = 40;

% taugrid is the old tau grid for quantile regression
% taugrid_midpoint is the mid point of tau segments
% taugrid_ue is the new uneven grid
[taugrid, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau);

% Part B) Load original data & keep WHITE observations
[X_orig, y_orig] = white_subsample(data);
clear data
% Calculate the number of observations in the subsamples
[nsample_orig, ncovar] = size(X_orig);
nsample = floor(1/n_subsample*nsample_orig);
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
recorder_qreg_repeat = nan(n_repeatsample,ncovar*ntau);
recorder_qreg_repeat_reshape = nan(n_repeatsample,ncovar,ntau);
recorder_WLS_repeat = nan(n_repeatsample,ncovar*ntau);

% piecewise constant result
recorder_WLS_start_repeat = nan(n_repeatsample,nvars);
fval_recorder_WLS_start_repeat = nan(1,n_repeatsample);
exit_recorder_WLS_start_repeat = nan(1,n_repeatsample);
recorder_WLS_start_repeat_sorted = nan(n_repeatsample,ncovar,ntau);

% piecewise linear result
recorder_pc_start_repeat = nan(n_repeatsample,nvars);
fval_recorder_pc_start_repeat = nan(1,n_repeatsample);
exit_recorder_pc_start_repeat = nan(1,n_repeatsample);
recorder_pc_start_repeat_recon = nan(n_repeatsample,ncovar,ntau);

% Mean of covariates
X_mean_repeat = nan(n_repeatsample,ncovar);

% Upper and lower bounds for parameters in the MLE estimation
lower=[zeros(1,(ncovar*ntau))-10000, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(ncovar*ntau))+10000), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
lower_pl=[zeros(1,(ncovar*ntau))-10, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper_pl=[(zeros(1,(ncovar*ntau))+100), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];

% Constants for MLE estimation, forcing the sum of the first two weights
% for mixtures not exceeding 1
b=1;
A = zeros(1,nvars);
A(ncovar*ntau+1) = 1;
A(ncovar*ntau+2) = 1;
    
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
    
    X_mean = mean(X);
    X_mean_repeat(j,:) = X_mean;
    
    % Part E) qreg and WLS
    [fit] = quantlsfVector(X,y,taugrid);
    fit_1 = reshape(fit, ncovar*ntau,1);
    recorder_qreg = fit_1';
    recorder_qreg_repeat(j,:) = recorder_qreg;
    recorder_qreg_repeat_reshape(j,:,:) = fit';
    recorder_WLS_repeat(j,:) =  WLS_step(fit,X,y,n_WLS_iter);
    
    % Part F) MLE
    start = [recorder_WLS_repeat(j,:), para_dist_default];
    [fit_hat,fval,exitflag] = fmincon(@(x)gradllfCovarparavector(x, ntau, nsample, nmixtures,1,y, X),start,A,b,[],[],lower,upper,[],options);   
    recorder_WLS_start_repeat(j,:) = fit_hat;
    fval_recorder_WLS_start_repeat(j) = fval;
    exit_recorder_WLS_start_repeat(j) = exitflag;
    
    
    % G) Piecewise linear MLE
    % G-1) Sort piecewise constant result

    recorder_WLS_start_reshape = reshape(fit_hat(1:ntau*ncovar),ntau,ncovar);
    [beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,ncovar);
    beta_WLS_start_sorted = beta_WLS_start_sorted';
    recorder_WLS_start_repeat_sorted(j,:,:) = beta_WLS_start_sorted;

    % G-2) Construct the start value
    [fit_1_temp] = construct_pl_start(beta_WLS_start_sorted, ncovar, ntau);
    start = [fit_1_temp, fit_hat([(end - 3*nmixtures + 3) : end])];
    
    [fit_hat,fval,exitflag,output] = fmincon(@(x)gradl_CDF_Lei_GA_ue(x,  taugrid_ue, nmixtures, y', X'), start, A, b,[],[],lower_pl,upper_pl,[],options);
    recorder_pc_start_repeat(j,:) = fit_hat;
    fval_recorder_pc_start_repeat(j) = fval;
    exit_recorder_pc_start_repeat(j) = exitflag;
    
    beta_pc_s = reconstruct_beta( (reshape(fit_hat(1,[1:(ncovar*ntau)]), [ntau, ncovar]))'  );
    recorder_pc_start_repeat_recon(j,:,:) = beta_pc_s;
    
    if mod(j,50) == 0     
        clear X y v_sample v_sample_sorted
        save (sprintf('RData_pl_%s_sub%d_ntau%d_white%d_ue1',postfix,n_subsample,ntau,n_repeatsample))
    end

end

clear X y X_orig y_orig v_sample v_sample_sorted
save (sprintf('RData_pl_%s_sub%d_ntau%d_white%d_ue1',postfix,n_subsample,ntau,n_repeatsample))
