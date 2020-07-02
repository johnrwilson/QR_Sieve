% This version
% This codes gives a start value to piecewise linear MLE simulations.
% This version uses GA, instead of fmincon
% Haoyang Liu
% 12/09/2018
function [betas] = QR_sieve(X, y, ntau, n_WLS_iter, upper, lower, para_dist_default, A, b, do_mle, n_batches, n_epochs, learning_rate)

% Generate the grid of knots for tau

% taugrid is the old tau grid for quantile regression
% taugrid_midpoint is the mid point of tau segments
% taugrid_ue is the new uneven grid

% Question: Currently taugrid_midpoint is not used. Is there another
% formulation where we use the midpoints and that's why it's calculated?

% TODO: Figure out how to avoid calculating this inside the function

% Note to John: taugrid is for the piecewise constant MLE and the WLS since 
% we have to sort the results, and an equal grid ensures an invariant log 
% likelihood. taugrid_ue is for the piecewise linear MLE calculations since
% grid intervals won't matter.
[taugrid, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau);

% Part B) Preallocation of result variables
nsample = size(X,1);
ncovar = size(X,2);
%number of mixture components, using mixture of normals
% TODO: unfreeze this parameter
nmixtures=3;
% nvars is the total number of parameters
nvars = ncovar*ntau+3*nmixtures-2;
% Note: the 3 is because for each mixture we estimate (1) the mean, (2) the
% weight, and (3) the std dev. The -2 is because one mean and one weight
% are pinned down by the other ones.
    
X_mean = mean(X);

% Part E) qreg and WLS
[fit] = quantlsfVector(X,y,taugrid);
WLS =  WLS_step(fit,X,y,n_WLS_iter);

if do_mle
    
    % Part F) MLE
    start = [WLS, para_dist_default];
    % TODO: Unfreeze the minimizing function
    % TODO: Decide if we want to save fval and exitflag
    options=optimoptions(@fmincon,'GradObj','on');
    [fit_hat] = fmincon(@(x)gradllfCovarparavector(x, ntau, nmixtures,1,y, X),start,A,b,[],[],lower,upper,[],options);
    disp("finished first MLE sgd")
    [loss, ~] = gradllfCovarparavector(fit_hat, ntau, nmixtures,1,y, X);
%     disp(loss);
%     disp(fit_hat);
        
    
    % G) Piecewise linear MLE
    % G-1) Sort piecewise constant result
    WLS_start_reshape = reshape(fit_hat(1:ntau*ncovar),ntau,ncovar);
    % TODO: Determine if we want to same V_WLS_start
    [beta_WLS_start_sorted] = sortbeta_1(X_mean,WLS_start_reshape,ntau,ncovar);
    beta_WLS_start_sorted = beta_WLS_start_sorted';

    % G-2) Construct the start value
    [fit_1_temp] = construct_pl_start(beta_WLS_start_sorted, ncovar, ntau);
    start = [fit_1_temp, fit_hat([(end - 3*nmixtures + 3) : end])];

    % TODO: Again, allow the user to specify their own minimizer
    [loss, grad] = gradl_CDF_Lei_GA_ue(start, taugrid_ue, nmixtures,y', X');
    disp("Loss and gradient at initial value")
    disp(loss)
    disp(grad)
    f = @(x,y,X) gradl_CDF_Lei_GA_ue(x, taugrid_ue, nmixtures, y', X');
    verbose = true;
    [fit_hat] = sgd(f, start, y, X, n_batches, n_epochs, learning_rate, verbose);
    
    betas = reconstruct_beta((reshape(fit_hat(1,[1:(ncovar*ntau)]), [ntau, ncovar]))');
%     disp(fit_hat(ncovar*ntau+1:end))
    disp("finished second MLE sgd")
    [loss, ~] = gradl_CDF_Lei_GA_ue(fit_hat, taugrid_ue, nmixtures,y', X');
    disp(loss);
    
else
    
    betas = WLS;
    
end