% This the code testing the gradient function
% It was written by both Haoyang and Lei
% 12/2/2017

clear all

%% Part A) Set up input to the likelihood function
% A - 1) Parameters

% Tau-vector
    Ntau=33;
    taustep=1/Ntau;
    taugrid=[0:taustep:1-taustep];

    N = 1000;
    nocovar = 5;
    para = [1/2 1 2];
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

    c = para(:,1);
    beta = para(:,[2:end])';
    beta_v = beta(:);
    cbeta = [c', beta_v'];
    cbeta = cbeta(:)';
    % Very important: cbeta begins with the 5 c's, followed by all the
    % beta's 
    
%% Part B) Calculate the analytical gradient
    [lf, g] = gradl_Lei([cbeta, ParameterDist([1,2,4,5,7:9])], taugrid, nmixtures, y, x);

%% Part C) Calculate the empirical gradients
% C-1) A small disturbance to c and beta
    for any = [1:170]
        any
        epsilon = 0.001;
        para_temp = [cbeta, ParameterDist([1,2,4,5,7:9])];
        para_temp(any) = para_temp(any) + epsilon;
        [lf_ temp]  = gradl_Lei(para_temp, taugrid, nmixtures, y, x);
        g_any(1,any) = (lf_-lf)/epsilon;
    end

    g_compare = [g(1:170); g_any];


save test_4_100

