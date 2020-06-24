% This function converts the qreg results to start values of MLE for
% piecewise linear
% Haoyang Liu
% 11/17/2017
function [c_start,beta_start] = pc_pl(beta_constants,ncovar,ntau)


beta_c_matrix = reshape(beta_constants,[ntau+1, ncovar])';

beta_l_matrix = (beta_c_matrix(:,[2:end]) - beta_c_matrix(:,[1:end-1]))*ntau;
    
beta_l_matrix = beta_l_matrix';

beta_start = beta_l_matrix(:)';


c_start = beta_c_matrix(:,1)';





end
