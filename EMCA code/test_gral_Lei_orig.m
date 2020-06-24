clear all

% Tau-vector
Ntau=33;
taustep=1/Ntau;
taugrid=[0:taustep:1-taustep];

%% Part A) Set up input to the likelihood function
% A - 1) Parameters
N = 100;
nocovar = 5;
para = [1/2 1 2];
para = repmat(para, 1, Ntau/3);
para = [1 para]; 

para = repmat(para, nocovar, 1); % each covariate has a set of beta values
para(2,:) = para(2,:)+1;
para(3,:) = para(3,:)-1;
%para(4,:) = para(4,:)+2;
%para(5,:) = para(5,:)-2;

nmixtures = 3;
ParameterDist = [[1/2 1/4 1/4],[-3 2 4],[1 1 1]];

% A - 2) Simulate x and y
y = rand(1, N);
x = rand(nocovar, N);

c = para(:,1);
beta = para(:,[2:end])';
beta_v = beta(:);
cbeta = [c', beta_v'];
cbeta = cbeta(:)';
[lf, g] = gradl_Lei_orig([cbeta, ParameterDist([1,2,4,5,7:9])], taugrid, nmixtures, y, x);


% C-3) A small disturbance to lambda
for any = [1:177]
epsilon = 0.001;
para_temp = [cbeta, ParameterDist([1,2,4,5,7:9])];
para_temp(any) = para_temp(any) + epsilon;
[lf_any, g] = gradl_Lei(para_temp, taugrid, nmixtures, y, x);
g_any(any) = (lf_any-lf)/epsilon;
end
g_compare = [g; g_any];
save test_4

