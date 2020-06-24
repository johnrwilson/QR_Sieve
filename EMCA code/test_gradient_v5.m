
clear all

% Tau-vector
Ntau=33;
taustep=1/Ntau;
taugrid=[0:taustep:1-taustep];

%% Part A) Set up input to the likelihood function
% A - 1) Parameters
N = 10000;
nocovar = 5;
para = [0 1/2 1];
para = repmat(para, 1, Ntau/3);
para = [1 para]; 

para = repmat(para, nocovar, 1); % each covariate has a set of beta values
para(2,:) = para(2,:)+1;
para(3,:) = para(3,:)-1;

nmixtures = 3;
ParameterDist = [[1/2 1/4 1/4],[-3 2 4],[1 1 1]];
lambda = ParameterDist(1:nmixtures);
mu = ParameterDist(nmixtures+1: 2*nmixtures);
sigma = ParameterDist(2*nmixtures+1: 3*nmixtures);

% A - 2) Simulate x and y
y = rand(1, N);
x = rand(nocovar, N);

%% Part B) Empirical Gradients for beta and c
leachsample = integral(@(t)fun(para, taugrid, t, taustep, x, y, ParameterDist, nmixtures), 0, 1, 'ArrayValued', true);
lf = sum(log(leachsample))/N;

g_empirical = nan(nocovar, (Ntau+1));

for j_covar = [1:nocovar] 
    for j_ntau = [1:(Ntau+1)]
        epsilon = 0.0001;
        para_temp = para;
        para_temp(j_covar,j_ntau) = para_temp(j_covar,j_ntau) + epsilon;
        leachsample_temp = integral(@(t)fun(para_temp, taugrid, t, taustep, x, y, ParameterDist, nmixtures), 0, 1,'ArrayValued', true)';
        lf_temp = sum(log(leachsample_temp))/N;
        g_empirical(j_covar,j_ntau) = (lf_temp-lf)/epsilon;
    end
    
end



%% Part C) Analytical gradients for c, beta, lambda, mu, sigma
phi_u = zeros(N, nmixtures);
u_phi_u = zeros(N, nmixtures);
u2_phi_u = zeros(N, nmixtures);

for j = 1:nmixtures
    phi_u(:, j)=integral(@(t)fun_phi_u(para, taugrid, t, taustep, x, y, mu(j), sigma(j)), 0, 1,'ArrayValued', true);
    u_phi_u(:, j)=integral(@(t)fun_u_phi_u(para, taugrid, t, taustep, x, y, mu(j), sigma(j)), 0, 1,'ArrayValued', true);
    u2_phi_u(:, j)=integral(@(t)fun_u2_phi_u(para, taugrid, t, taustep, x, y, mu(j), sigma(j)), 0, 1,'ArrayValued', true);
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
    newtau = tau_lower:0.001:tau_lower+taustep;
    newtaugrid = tau_lower:taustep:tau_lower+taustep;
    newtaugrid2 = tau_lower:taustep:taugrid(Ntau);
    k(:,:,i) = Ntau-i;
    
    for j = 1:nocovar
    gbeta1eachsample(j,:,i) = integral(@(t)fun_beta_1(para, newtaugrid, t, taustep, ...
                  x, x(j,:), y, tau_lower, ParameterDist, nmixtures), tau_lower, tau_lower+taustep,'ArrayValued', true)';
    gbeta2eachsample(j,:,i) = integral(@(t)fun_beta_2(para, newtaugrid2, t, taustep, ...
                  x, x(j,:), y, ParameterDist, nmixtures), tau_lower+taustep, 1,'ArrayValued', true)';
    end
  
end
gbetaeachsample = gbeta1eachsample + gbeta2eachsample;
gbetaeachsample = gbetaeachsample./permute(repmat(leachsample, [Ntau,1, 5]),[3, 2, 1]);

% gradients for c1
gceachsample = zeros(nocovar, N);
for j = 1:nocovar
 gceachsample(j,:) = integral(@(t)fun_c(para, taugrid, t, taustep, ...
                  x, x(j,:), y, ParameterDist, nmixtures), 0, 1,'ArrayValued', true, 'AbsTol', 0.001)';
end
gceachsample = gceachsample./repmat(leachsample, nocovar, 1);


glambda = sum(glambdaeachsample)/N;
gmu = sum(gmueachsample)/N;
gsigma = sum(gsigmaeachsample)/N;
gbeta = squeeze(sum(gbetaeachsample,2))/N;
gc = sum(gceachsample, 2)/N;


g_compare = [g_empirical; gc, gbeta];
% b(t)
function value = beta(theta,tgrid,t,dtau)
    t = t';
    tgrid = repmat(tgrid, length(t), 1);
% 	tgrid is the grid of taus
% 	tvector is the vector of taus that we want to evaluate
% 	theta is the vector of parameters
% 	dtau is step size
	[tbelow,tbelow_index]=max((tgrid<t).*tgrid, [], 2);
    thetamatrix = repmat(theta, length(t), 1);
    indexmat = ones(length(t), length(theta));
    indexmat = (cumsum(indexmat,2)<=tbelow_index);
    indexmat(:,1) = 0 ;
    constantterm = sum(thetamatrix.*indexmat,2)*dtau;
	value = theta(1)+constantterm+theta(tbelow_index+1)'.*(t-tbelow);
end

% f(y-xb(t))
function value = ftau(r_matrix, ParameterDist, nmixtures);
    lambda = ParameterDist(1:nmixtures)';
    mu = ParameterDist(nmixtures+1: 2*nmixtures)';
    sigma = ParameterDist(2*nmixtures+1: 3*nmixtures)';
    
    z = 0*r_matrix;
    for j = 1:nmixtures
        z = z+lambda(j)/sigma(j)*normpdf((r_matrix-mu(j))/sigma(j));
    end
    value = z;
end

% erroreachsample
function value = error(theta, tgrid, t, dtau, x, y)
    value = y;
    k = 1;
    while k <= size(x, 1)
        value = value-x(k,:)*beta(theta(k,:), tgrid, t, dtau);
        k = k+1;
    end
end

% function used to calcualte the integral of f(y-xb(t))
function value = fun(theta, tgrid, t, dtau, x, y, ParameterDist, nmixtures)
    erroreachsample = error(theta, tgrid, t, dtau, x, y);
    value = ftau(erroreachsample, ParameterDist, nmixtures);
end

function value = fun_phi_u(theta, tgrid, t, dtau, x, y, mu, sigma)
    erroreachsample = error(theta, tgrid, t, dtau, x, y);
    value = normpdf((erroreachsample-mu)/sigma);
end

function value = fun_u_phi_u(theta, tgrid, t, dtau, x, y, mu, sigma)
    erroreachsample = error(theta, tgrid, t, dtau, x, y);
    value = ((erroreachsample-mu)/sigma).*normpdf((erroreachsample-mu)/sigma);
end

function value = fun_u2_phi_u(theta, tgrid, t, dtau, x, y, mu, sigma)
    erroreachsample = error(theta, tgrid, t, dtau, x, y);
    value = (((erroreachsample-mu)/sigma).^2).* normpdf((erroreachsample-mu)/sigma);
end

% calculate the sum value inside beta integral
function value = fsum(r_matrix, ParameterDist, nmixtures);
    lambda = ParameterDist(1:nmixtures)';
    mu = ParameterDist(nmixtures+1: 2*nmixtures)';
    sigma = ParameterDist(2*nmixtures+1: 3*nmixtures)';
    
    z = 0*r_matrix;
    for j = 1:nmixtures
        z = z+lambda(j)/(sigma(j).^2).*((r_matrix-mu(j))/sigma(j)).*normpdf((r_matrix-mu(j))/sigma(j));
    end
    value = z;
end

function value = fun_beta_1(theta, tgrid, t, dtau, x, x2, y, taulower, ParameterDist, nmixtures)
    erroreachsample = error(theta, tgrid, t, dtau, x, y);
    value = fsum(erroreachsample, ParameterDist, nmixtures).*x2.*(t-taulower);
end

function value = fun_beta_2(theta, tgrid, t, dtau, x, x2, y, ParameterDist, nmixtures)
    erroreachsample = error(theta, tgrid, t, dtau, x, y);
    value = fsum(erroreachsample, ParameterDist, nmixtures).*x2.*dtau;
end

function value = fun_c(theta, tgrid, t, dtau, x, x2, y, ParameterDist, nmixtures)
    erroreachsample = error(theta, tgrid, t, dtau, x, y);
    value = fsum(erroreachsample, ParameterDist, nmixtures).*x2;
end




