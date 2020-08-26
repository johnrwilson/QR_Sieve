load('../../results_seeded.mat', 'taugrid_ue', 'betas_mle');
load('results_sgd_seeded.mat', 'betas_sgd');
load('results_1mle_seeded.mat', 'betas_1mle');

% betas_1mle = reshape(

figure()
hold on;
title("Beta 1")
plot(taugrid_ue, sort(betas_mle(1,:)), 'LineWidth', 2);
plot(taugrid_ue, sort(betas_1mle(1,:)), 'LineWidth', 2);
plot(taugrid_ue, sort(betas_sgd(1,:)), 'LineWidth', 2);
legend("GA", "Piecewise constant only", "SGD");
% legend("GA", "SGD");



figure()
hold on;
title("Beta 2")
plot(taugrid_ue, sort(betas_mle(2,:)), 'LineWidth', 2);
plot(taugrid_ue, sort(betas_1mle(2,:)), 'LineWidth', 2);
plot(taugrid_ue, sort(betas_sgd(2,:)), 'LineWidth', 2);
legend("GA", "Piecewise constant only", "SGD");
% legend("GA", "SGD");

figure()
hold on;
title("Beta 3")
plot(taugrid_ue, sort(betas_mle(3,:)), 'LineWidth', 2);
plot(taugrid_ue, sort(betas_1mle(3,:)), 'LineWidth', 2);
plot(taugrid_ue, sort(betas_sgd(3,:)), 'LineWidth', 2);
legend("GA", "Piecewise constant only", "SGD");
% legend("GA", "SGD");