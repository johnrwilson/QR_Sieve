function plot_bootstrap(betas, fit_hat, fit_bootstrap, betas_bootstrap)

[ncovar, ntau, bootstrap] = size(betas_bootstrap);

betas_std = std(betas_bootstrap, 0, 3);
[~, ~, taugrid_ue] = calculate_grid(ntau);
disp(ntau)

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

nmixtures = (size(fit_bootstrap,2) + 2 - ncovar * ntau)/3;

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

density_dist = zeros(bootstrap,n_points);

parfor j = 1:bootstrap
    lambdas_short = fit_bootstrap(j,ncovar*ntau+1:ncovar*ntau+nmixtures-1).^2;
    lambdas = [lambdas_short, 1 - sum(lambdas_short)]
    mus_short = fit_bootstrap(j,ncovar*ntau+nmixtures:ncovar*ntau+2*nmixtures-2);
    mus = [mus_short, -sum(lambdas_short.*mus_short)/lambdas(end)]
    sigmas = fit_bootstrap(j,end-2:end)
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