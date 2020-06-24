% This code runs quantile regression and WLS on the GDP growth data  

% 12/17/2018
clear
close all 



data = csvread('gdp_nfci.csv',1,0);

gdp_growth = (((data([5:end],1) - data([1:end-4],1))./data([1:end-4],1))*100);
gdp_growth_future = gdp_growth([5:end]);

nfci = data([1:end-8],2)';
nsample = length(nfci);
ncovar = 2;
ntau = 19;
[taugrid, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau);

 X = [ones(1,nsample); gdp_growth(1:end-4)'; nfci]';
  X = [ gdp_growth(1:end-4)'; nfci]';
[fit_qreg] = quantlsfVector(X,gdp_growth_future,taugrid);
n_WLS_iter = 50;
[temp,fit_WLS_1] =  WLS_step(fit_qreg,X,gdp_growth_future,n_WLS_iter);
n_WLS_iter = 200;
[temp,fit_WLS_2] =   WLS_step(fit_qreg,X,gdp_growth_future,n_WLS_iter);
n_WLS_iter = 500;
[temp,fit_WLS_3] =   WLS_step(fit_qreg,X,gdp_growth_future,n_WLS_iter);

fit_qreg_new = quantlsfVector(X,gdp_growth_future,taugrid_midpoint);

fit_compare = [fit_qreg';fit_WLS_1';fit_WLS_2';fit_WLS_3';fit_qreg_new'];


nvars = ncovar*ntau + 7;
b=1;
A = zeros(1,nvars);
A(ncovar*ntau+1) = 1;
A(ncovar*ntau+2) = 1;
nmixtures = 3;
para_dist_default = [[1/3,1/3],[-1,0],[1,1,1]];


lower=[zeros(1,(ncovar*ntau))-20, (zeros(1,(nmixtures-1))+0.001),(zeros(1,(nmixtures-1))-10.01),(zeros(1,nmixtures)+0.01)];
upper=[(zeros(1,(ncovar*ntau))+20), (zeros(1,(nmixtures-1))*0+1),(zeros(1,(nmixtures-1))+10),ones(1,nmixtures)*10];
 options=optimoptions(@fmincon,'GradObj','on');
start = [temp, para_dist_default];
[fit_hat,fval,exitflag] = fmincon(@(x)gradllfCovarparavector(x, ntau, nsample, nmixtures,1,gdp_growth_future, X),start,A,b,[],[],lower,upper,[],options);


fit_compare = [fit_qreg';fit_WLS_1';fit_WLS_2';fit_WLS_3';fit_qreg_new';reshape(fit_hat(1:ntau*ncovar),ntau,ncovar)'];




