clear all

% Tau-vector
Ntau=33;
taustep=1/Ntau;
taugrid=[0:taustep:1-taustep];

%% Part A) Set up input to the likelihood function
% A - 1) Parameters
N = 100;
nocovar = 5;
para = [0 1/2 1];
para = repmat(para, 1, Ntau/3);
para = [1 para]; 

para = repmat(para, nocovar, 1); % each covariate has a set of beta values
para(2,:) = para(2,:)+1;
para(3,:) = para(3,:)-1;

nmixtures = 3;
ParameterDist = [[1/2 1/4 1/4],[-3 2 4],[1 1 1]];

% A - 2) Simulate x and y
y = rand(1, N);
x = rand(nocovar, N);


[lf, g_temp] = gradl_nopre_Lei(para, ParameterDist, taugrid, nmixtures, nocovar, N, y, x);

%% Part B) Test the gradient function
% C-3) A small disturbance to lambda
for any = [1:177]
epsilon = 0.001;
para_temp = para;
para_temp(any) = para_temp(any) + epsilon;
[lf_any, g] = gradl_Lei(para_temp, taugrid, nmixtures, y, x);
g_any(any) = (lf_any-lf)/epsilon;
end
g_compare = [g; g_any];



return;
% C-3) A small disturbance to lambda
epsilon = 0.001;
para_temp = ParameterDist;
para_temp(1,1) = para_temp(1,1) + epsilon;
[lf_lambda, g_temp] = gradl_nopre_Lei(para, para_temp, taugrid, nmixtures, nocovar, N, y, x);
g_lambda = (lf_lambda-lf)/epsilon;

% C-4) A small disturbance to mu
epsilon = 0.001;
para_temp = ParameterDist;
para_temp(1,4) = para_temp(1,4) + epsilon;
[lf_mu, g_temp] = gradl_nopre_Lei(para, para_temp, taugrid, nmixtures, nocovar, N, y, x);
g_mu = (lf_mu-lf)/epsilon;

% C-4) A small disturbance to sigma
epsilon = 0.001;
para_temp = ParameterDist;
para_temp(1,7) = para_temp(1,7) + epsilon;
[lf_sigma, g_temp] = gradl_nopre_Lei(para, para_temp, taugrid, nmixtures, nocovar, N, y, x);
g_sigma = (lf_sigma-lf)/epsilon;
