% This function converts the MLE slope results to piecewise linear level
% Haoyang Liu
% 1/23/2018
function [beta_pc] = pl_pc(beta_constants, tau_step)

ntau = 1/tau_step;
beta_pc(:,1) = beta_constants(:,1);

for j_tau = [1 : (ntau)+0.0001]
    beta_pc(:,j_tau+1) = beta_pc(:,j_tau) + tau_step*beta_constants(:,j_tau+1);
end


