function [betas, fit_hat, betas_bootstrap, fit_hat_bootstrap] = ...
    sieve_mle(X, y, bootstrap, ntau, nmixtures, n_WLS_iter, lower, upper, ...
    optimizer, make_plot)

if ~exist('make_plot','var') || isempty(make_plot)
     % last parameter does not exist, so default it to something
     make_plot = false;
end

if ~exist('optimizer','var') || isempty(optimizer)
     % last parameter does not exist, so default it to something
     n_batches = 30;
     n_epochs = 1000;
     learning_rate = .00001;
     decay = .999;
     verbose = true;
     optimizer = {'SGD', n_batches, n_epochs, learning_rate, decay, verbose};
end

if ~exist('do_mle','var') || isempty(do_mle)
     % last parameter does not exist, so default it to something
     do_mle = true;
end

std_y = std(y);
if ~exist('upper','var') || isempty(upper)
     % last parameter does not exist, so default it to something
     E_y = mean(y);
     upper = [10000, 1, E_y + 3 * std_y, 10*std_y]
end

if ~exist('lower','var') || isempty(lower)
    % last parameter does not exist, so default it to something
    E_y = mean(y);
    lower = [-10000, .0001, -E_y - 3 * std_y, .01*std_y]
end

if ~exist('n_WLS_iter','var') || isempty(n_WLS_iter)
     % last parameter does not exist, so default it to something
     n_WLS_iter = 40;
end

if ~exist('nmixtures','var') || isempty(nmixtures)
     % last parameter does not exist, so default it to something
     nmixtures = 3;
end

if ~exist('ntau','var') || isempty(ntau)
     % last parameter does not exist, so default it to something
     ntau = 15;
end

% Generate the grid of knots for tau

% taugrid is the old tau grid for quantile regression
% taugrid_ue is the new uneven grid

[taugrid, ~, taugrid_ue] = calculate_grid(ntau);

nsample = size(X,1);
ncovar = size(X,2);

% nvars is the total number of parameters
nvars = ncovar*ntau+3*nmixtures-2;
% Note: the 3 is because for each mixture we estimate (1) the mean, (2) the
% weight, and (3) the std dev. The -2 is because one mean and one weight
% are pinned down by the other ones.

% Constants
b=1;

A = zeros(1,nvars);
A(ncovar*ntau+1:ncovar*ntau+nmixtures-1) = 1;

lower=[zeros(1,(ncovar*ntau))+lower(1), (zeros(1,(nmixtures-1))+lower(2)), ...
    (zeros(1,(nmixtures-1))+lower(3)),(zeros(1,nmixtures)+lower(4))];

upper=[zeros(1,(ncovar*ntau))+upper(1), (zeros(1,(nmixtures-1))+upper(2)), ...
    (zeros(1,(nmixtures-1))+upper(3)),(zeros(1,nmixtures)+upper(4))];

% define the starting values for the guess of distributional parameters
if nmixtures == 1

elseif mod(nmixtures, 2) == 0
    lambda_start = repmat([1/nmixtures], 1, nmixtures-1);
    sigma_start = repmat([sqrt(.75) * std_y], 1, nmixtures);
    mu_start = [repmat([-1], 1, nmixtures/2), ...
        repmat([1], 1, nmixtures/2-1)] * sqrt(.75) * std_y;
    para_dist_default = [lambda_start, mu_start, sigma_start];
else
    lambda_start = repmat([1/nmixtures], 1, nmixtures-1);
    sigma_start = repmat([sqrt(.75) * std_y], 1, nmixtures);
    mu_start = [repmat([-1], 1, (nmixtures-1)/2), 0, ...
        repmat([1], 1, (nmixtures-1)/2-1)] * sqrt(.75) * std_y;
    para_dist_default = [lambda_start, mu_start, sigma_start];
end

X_mean = mean(X);

% Part E) qreg and WLS
[fit] = quantlsfVector(X,y,taugrid);
fit = reshape(fit, ncovar*ntau, 1)';
% WLS =  WLS_step(fit,X,y,n_WLS_iter);

% Part F) MLE
start = [fit, para_dist_default];
% TODO: Unfreeze the minimizing function
%     n_batches = 50;
%     n_epochs = 50;
%     learning_rate = .00001;
%     decay = .999;
%     verbose = true;

%      f = @(x,y,X)gradllfCovarparavector(x, ntau, nmixtures,1,y, X);
%     [fit_hat] = sgd(f, start, y, X, n_batches, n_epochs, learning_rate, decay, verbose);

options=optimoptions(@fmincon,'GradObj','on');
[fit_hat] = fmincon(@(x)gradllfCovarparavector(x, ntau, nmixtures,1,y, X),...
    start,A,b,[],[],lower,upper,[],options);
disp(gradllfCovarparavector(fit_hat, ntau, nmixtures,1,y, X));
disp("Finished first MLE")

% G) Piecewise linear MLE
% G-1) Sort piecewise constant result
WLS_start_reshape = reshape(fit_hat(1:ntau*ncovar),ntau,ncovar);
[beta_WLS_start_sorted] = sortbeta_1(X_mean,WLS_start_reshape,ntau,ncovar);
beta_WLS_start_sorted = beta_WLS_start_sorted';

% G-2) Construct the start value
[fit_1_temp] = construct_pl_start(beta_WLS_start_sorted, ncovar, ntau);
lambdas = sqrt(fit_hat([(end - 3*nmixtures + 3) : (end - 2*nmixtures + 1)]));
start = [fit_1_temp, lambdas, fit_hat([(end - 2*nmixtures + 2) : end])];

if optimizer{1} == "GA"

    disp("Solving second MLE with genetic algorithm")

    opts = optimizer{2};
    opts.InitialPopulationMatrix = start;
    [fit_hat] = ga(@(x)gradl_CDF_Lei_GA_ue(x, taugrid_ue, nmixtures, y', X'), nvars, A, b,[],[],lower,upper,[],opts);

elseif optimizer{1} == "SGD"

    disp("Solving second MLE with stochastic gradient descent")

    n_batches = optimizer{2};
    n_epochs = optimizer{3};
    learning_rate = optimizer{4};
    decay = optimizer{5};
    verbose = optimizer{6};

    f = @(x,y,X) gradl_CDF_Lei_GA_ue_free_lambdas(x, taugrid_ue, nmixtures, y', X');
    [fit_hat] = sgd(f, start, y, X, n_batches, n_epochs, learning_rate, decay, verbose);

elseif optimizer{1} == "SA"

    disp("Solving second MLE with simulated annealing")

    %TODO: List requirements for custom optimizer here

    f = @(x) gradl_CDF_Lei_GA_ue(x, taugrid_ue, nmixtures, y', X');
    fit_hat = simulannealbnd(f, start, lower, upper, optimizer{2});

elseif optimizer{1} == "FM"

    disp("Solving second MLE with fminsearch")

    %TODO: List requirements for custom optimizer here

    f = @(x) gradl_CDF_Lei_GA_ue(x, taugrid_ue, nmixtures, y', X');
    fit_hat = fminsearch(f, start, optimizer{2});

end

betas = reconstruct_beta((reshape(fit_hat(1,[1:(ncovar*ntau)]), [ntau, ncovar]))');
disp("Finished second MLE")

fit_hat_bootstrap = fit_hat;
betas_bootstrap = betas;
    
if bootstrap > 1
    if optimizer{1} == "SGD"
        optimizer{6} = false;
    end
    fit_hat_bootstrap = zeros(bootstrap, nvars);
    betas_bootstrap = zeros(ncovar, ntau, bootstrap);
    parfor i = 1:bootstrap
        sample_index = randsample(1:nsample, nsample, true);
        y_sample = y(sample_index);
        X_sample = X(sample_index, :);
        [betas_subsample, fit_hat_subsample] = QR_sieve(X, y, 1, ntau, ...
            nmixtures, n_WLS_iter, lower, upper, optimizer, make_plot);
        fit_hat_bootstrap(i,:) = fit_hat_subsample
        betas_bootstrap(:,:,i) = betas_subsample
    end
    if make_plot
        betas_std = std(betas_bootstrap, 0, 3);
        for i = 1:ncovar
            figure()
            hold on
            plot(taugrid_ue, betas(i,:), 'LineWidth', 2)
            plot(taugrid_ue, betas(i,:) + 1.96 * betas_std(i,:), 'LineWidth', 2)
            plot(taugrid_ue, betas(i,:) - 1.96 * betas_std(i,:), 'LineWidth', 2)
            xlabel('$\tau$', 'interpreter', 'latex', 'FontSize', 16)
            ylabel_text = sprintf("$\\beta_{%d}(\\tau)$", i);
            ylabel(ylabel_text, 'interpreter', 'latex', 'FontSize', 16);
        end
        
        lambdas_short = fit_hat(ncovar*ntau+1:ncovar*ntau+nmixtures-1).^2;
        lambdas = [lambdas_short, 1 - sum(lambdas_short)];
        mus_short = fit_hat(ncovar*ntau+nmixtures:ncovar*ntau+2*nmixtures-2);
        mus = [mus_short, -sum(lambdas_short.*mus_short)/lambdas(end)];
        sigmas = fit_hat(nvars-2:nvars);
        
        dom_min = min(mus - 3 * sigmas);
        dom_max = max(mus + 3 * sigmas);
        n_points = 100;
        dom = linspace(dom_min, dom_max, n_points);
        
        density_y = zeros(1,n_points);
        
        for i = 1:nmixtures
            density_y = density_y + lambdas(i) * normpdf(dom, mus(i), sigmas(i));
        end
        
        density_dist = zeros(bootstrap,n_points);
        
        parfor j = 1:bootstrap
            lambdas_short = fit_hat_bootstrap(j,ncovar*ntau+1:ncovar*ntau+nmixtures-1).^2;
            lambdas = [lambdas_short, 1 - sum(lambdas_short)];
            mus_short = fit_hat_bootstrap(j,ncovar*ntau+nmixtures:ncovar*ntau+2*nmixtures-2);
            mus = [mus_short, -sum(lambdas_short.*mus_short)/lambdas(end)]
            sigmas = fit_hat_bootstrap(j,nvars-2:nvars)
            for i = 1:nmixtures
                density_dist(j,:) = density_dist(j,:) + lambdas(i) * normpdf(dom, mus(i), sigmas(i));
            end
        end
        
        density_std = std(density_dist);
        
        figure()
        hold on;
        plot(dom, density_y)
        xlabel('Measurement Error', 'FontSize', 16)
        ylabel('Density', 'FontSize', 16);
        plot(dom, density_y + 1.96 * density_std, 'LineWidth', 2)
        plot(dom, density_y - 1.96 * density_std, 'LineWidth', 2)
    end
    
else
    if make_plot
        for i = 1:ncovar
            figure()
            plot(taugrid_ue, betas(i,:), 'LineWidth', 2)
            xlabel('$\tau$', 'interpreter', 'latex', 'FontSize', 16)
            ylabel_text = sprintf("$\\beta_{%d}(\\tau)$", i);
            ylabel(ylabel_text, 'interpreter', 'latex', 'FontSize', 16);
        end
        
        lambdas_short = fit_hat(ncovar*ntau+1:ncovar*ntau+nmixtures-1).^2;
        lambdas = [lambdas_short, 1 - sum(lambdas_short)];
        mus_short = fit_hat(ncovar*ntau+nmixtures:ncovar*ntau+2*nmixtures-2);
        mus = [mus_short, -sum(lambdas_short.*mus_short)/lambdas(end)];
        sigmas = fit_hat(end-2:end);
        
        dom_min = min(mus - 3 * sigmas);
        dom_max = max(mus + 3 * sigmas);
        n_points = 100;
        dom = linspace(dom_min, dom_max, n_points);
        
        density_y = zeros(1,n_points);
        
        for i = 1:nmixtures
            density_y = density_y + lambdas(i) * normpdf(dom, mus(i), sigmas(i));
        end
        
        figure()
        plot(dom, density_y, 'LineWidth', 2)
        xlabel('Measurement Error', 'FontSize', 16)
        ylabel('Density', 'FontSize', 16);
    end
end