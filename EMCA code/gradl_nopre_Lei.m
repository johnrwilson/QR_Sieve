% This function implements the gradient function by Lei
% para is a matrix. The first column is the c vector, and the rest of the
% matrix are betas
function [lf, g] = gradl_nopre_Lei(para, ParameterDist, taugrid, nmixtures, nocovar, N, yt, X)

% Tau step. Be careful. If the tau grid gets non-uniform, I need to ask Lei to change the function. 
taustep = taugrid(2) - taugrid(1);
Ntau = length(taugrid);


leachsample = integral(@(t)fun(para, taugrid, t, taustep, X, yt, ParameterDist, nmixtures), 0, 1, 'ArrayValued', true);

lf = sum(log(leachsample))/N;

lambda = ParameterDist(1:nmixtures);
mu = ParameterDist(nmixtures+1: 2*nmixtures);
sigma = ParameterDist(2*nmixtures+1: 3*nmixtures);

phi_u = zeros(N, nmixtures);
u_phi_u = zeros(N, nmixtures);
u2_phi_u = zeros(N, nmixtures);

for j = 1:nmixtures
    phi_u(:, j)=integral(@(t)fun_phi_u(para, taugrid, t, taustep, X, yt, mu(j), sigma(j)), 0, 1,'ArrayValued', true);
    u_phi_u(:, j)=integral(@(t)fun_u_phi_u(para, taugrid, t, taustep, X, yt, mu(j), sigma(j)), 0, 1,'ArrayValued', true);
    u2_phi_u(:, j)=integral(@(t)fun_u2_phi_u(para, taugrid, t, taustep, X, yt, mu(j), sigma(j)), 0, 1,'ArrayValued', true);
end

glambdaeachsample = phi_u./repmat(sigma, N, 1)./repmat(leachsample', 1, nmixtures);
gmueachsample = u_phi_u.*repmat(lambda, N,1)./repmat(sigma.^2, N,1)./repmat(leachsample', 1, nmixtures);
gsigma1eachsample = -1*phi_u.*repmat(lambda, N,1)./repmat(sigma.^2, N,1)./repmat(leachsample', 1, nmixtures);
gsigma2eachsample = u2_phi_u.*repmat(lambda, N,1)./repmat(sigma.^2, N,1)./repmat(leachsample', 1, nmixtures);
gsigmaeachsample = gsigma1eachsample + gsigma2eachsample;

% gradients for beta
gbeta1eachsample = zeros(nocovar, N, Ntau);
gbeta2eachsample = zeros(nocovar, N, Ntau);

for i = 1:Ntau
    tau_lower = taugrid(i);
    newtaugrid = tau_lower:taustep:tau_lower+taustep;
    newtaugrid2 = tau_lower:taustep:taugrid(Ntau);
    
    for j = 1:nocovar

    gbeta1eachsample(j,:,i) = integral(@(t)fun_beta_1(para, newtaugrid, t, taustep, ...
                  X, X(j,:), yt, tau_lower, ParameterDist, nmixtures), tau_lower, tau_lower+taustep,'ArrayValued', true)';     
    gbeta2eachsample(j,:,i) = integral(@(t)fun_beta_2(para, newtaugrid2, t, taustep, ...
                  X, X(j,:), yt, ParameterDist, nmixtures), tau_lower+taustep, 1,'ArrayValued', true)';
    end
end

gbetaeachsample = gbeta1eachsample + gbeta2eachsample;
gbetaeachsample = gbetaeachsample./permute(repmat(leachsample, [Ntau,1, nocovar]),[3, 2, 1]);


% gradients for c1
gceachsample = zeros(nocovar, N);
for j = 1:nocovar
 gceachsample(j,:) = integral(@(t)fun_c(para, taugrid, t, taustep, ...
                  X, X(j,:), yt, ParameterDist, nmixtures), 0, 1,'ArrayValued', true, 'AbsTol', 0.001)';
end
gceachsample = gceachsample./repmat(leachsample, nocovar, 1);


glambda = sum(glambdaeachsample)/N;
gmu = sum(gmueachsample)/N;
gsigma = sum(gsigmaeachsample)/N;
gbeta = squeeze(sum(gbetaeachsample,2))'/N;
gc = sum(gceachsample, 2)/N;


g_c_beta_v = [gc', (gbeta(:))'];


g = [g_c_beta_v, glambda, gmu, gsigma];













