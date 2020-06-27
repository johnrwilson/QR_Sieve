function [X_orig, y_orig] = white_subsample(data)

y_orig_1 = data(:,1);
X_orig_1 = [ones(size(y_orig_1)), data(:,2:end)];
[nsample_orig_1, ncovar_orig_1] = size(X_orig_1);

% Keep WHITE observations
y_orig = y_orig_1(X_orig_1(:,ncovar_orig_1)==0);
X_orig = X_orig_1(X_orig_1(:,ncovar_orig_1)==0,:);
X_orig = X_orig(:,[1:(end-1)]);
