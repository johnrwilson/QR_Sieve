X_mean = mean(X);
X_mean_repeat(j,:) = X_mean;

% Part E) qreg and WLS
[fit] = quantlsfVector(X,y,taugrid);
fit_1 = reshape(fit, ncovar*ntau,1);
recorder_qreg = fit_1';
recorder_qreg_repeat(j,:) = recorder_qreg;
recorder_qreg_repeat_reshape(j,:,:) = fit';
recorder_WLS_repeat(j,:) =  WLS_step(fit,X,y,n_WLS_iter);

% Part F) MLE
start = [recorder_qreg, para_dist_default];
[fit_hat,fval,exitflag] = fmincon(@(x)gradllfCovarparavector(x, ntau, nsample, nmixtures,1,y, X),start,A,b,[],[],lower,upper,[],options);
recorder_WLS_start_repeat(j,:) = fit_hat;
fval_recorder_WLS_start_repeat(j) = fval;
exit_recorder_WLS_start_repeat(j) = exitflag;


% G) Piecewise linear MLE
% G-1) Sort piecewise constant result

recorder_WLS_start_reshape = reshape(fit_hat(1:ntau*ncovar),ntau,ncovar);
[beta_WLS_start_sorted,V_WLS_start] = sortbeta_1(X_mean,recorder_WLS_start_reshape,ntau,ncovar);
beta_WLS_start_sorted = beta_WLS_start_sorted';
recorder_WLS_start_repeat_sorted(j,:,:) = beta_WLS_start_sorted;

% G-2) Construct the start value
[fit_1_temp] = construct_pl_start(beta_WLS_start_sorted, ncovar, ntau);
start = [fit_1_temp, fit_hat([(end - 3*nmixtures + 3) : end])];

opts.InitialPopulationMatrix = start;

[fit_hat,fval,exitflag,output] = ga(@(x)gradl_CDF_Lei_GA_ue(x,  taugrid_ue, nmixtures, y', X'), nvars, A, b,[],[],lower_pl,upper_pl,[],opts);
recorder_pc_start_repeat(j,:) = fit_hat;
fval_recorder_pc_start_repeat(j) = fval;
exit_recorder_pc_start_repeat(j) = exitflag;

beta_pc_s = reconstruct_beta( (reshape(fit_hat(1,[1:(ncovar*ntau)]), [ntau, ncovar]))'  );
recorder_pc_start_repeat_recon(j,:,:) = beta_pc_s;
