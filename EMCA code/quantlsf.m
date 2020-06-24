function b = quantlsf(X,y,q)
% Calculates Quantile Regression coefficients.
% Written by Moshe Buchinski and Gary Chamberlain.
% Alberto Abadie modified line 12.

limit=100;
[N, K] = size(X);

bls = X\y;
e = y-X*bls;
e = sort(e); 
bls(1) = bls(1) + e(round(N*q));
y = y-X*bls;

labelout = [1:K];
labelin = [K+1:K+N];
c = q*ones(1,N);

ind = find(y < 0);
if length(ind)>0,
y(ind) = -y(ind);
X(ind,:) = - X(ind,:);
labelin(ind) = -labelin(ind);
c(ind) = -c(ind) + 1.0;
end

%size(X)
T = [X, y; c*X, c*y];
clear X y c;
criterion = [T(N+1,1:K), -T(N+1,1:K)];
[margcost, s] = max(criterion);
count = 0;
for i = 1:K,
  if s > K,
     s = s-K; T(1:N,s) = -T(1:N,s);
     T(N+1,s) = margcost;
     labelout(s) = -labelout(s);
  end 
  ind = find(T(1:N,s) > 0 & abs(labelin') > K);
if length(ind)>0,
  [val, loc] = sort(T(ind,K+1) ./ T(ind,s));
clear val;
  r = ind(loc(1));
  init = 1;
  while (margcost - T(r,s)) > 0 & init < length(ind),
    T(N+1,:) = T(N+1,:) - T(r,:);
    margcost = T(N+1,s);
    T(r,:) = -T(r,:);
    labelin(r) = -labelin(r);
    init = init + 1;
    r = ind(loc(init));
  end
  T(:,K+2) = ([1:(N+1)]' == r);
  T(r,:) = T(r,:)/T(r,s);
  ind = find([1:N+1] ~= r);
  T(ind,:) = T(ind,:) - (T(ind,s)/T(r,s))*T(r,:);
  T(:,s) = T(:,K+2);
  move = labelout(s);
  labelout(s) = labelin(r);
  labelin(r) = move;
end

  if i < K,
    ind = find(abs(labelout) <= K);
    criterion = [T(N+1,ind), -T(N+1,ind)];
    [margcost, loc] = max(criterion);
    if loc <= length(ind),
      s = ind(loc);
    else
      s = K + ind(loc - length(ind));
    end
    count = count + 1;
  else
    criterion = [T(N+1,1:K), (-T(N+1,1:K) - 1.0)];
    [margcost, s] = max(criterion);
    count = count + 1;
  end
end
%
while margcost > 10*eps & count <= limit,
  if s > K,
     s = s-K; T(1:N,s) = -T(1:N,s);
     T(N+1,s) = margcost;
     labelout(s) = -labelout(s);
  end 
  ind = find(T(1:N,s) > 0 & abs(labelin') > K);
  [val, loc] = sort(T(ind,K+1) ./ T(ind,s));
clear val;
  r = ind(loc(1));
  init = 1;
  while (margcost - T(r,s)) > 0 & init < length(ind),
    T(N+1,:) = T(N+1,:) - T(r,:);
    margcost = T(N+1,s);
    T(r,:) = -T(r,:);
    labelin(r) = -labelin(r);
    init = init + 1;
    r = ind(loc(init));
  end
  T(:,K+2) = ([1:(N+1)]' == r);
  T(r,:) = T(r,:)/T(r,s);
  ind = find([1:N+1] ~= r);
  T(ind,:) = T(ind,:) - (T(ind,s)/T(r,s))*T(r,:);
  T(:,s) = T(:,K+2);
  move = labelout(s);
  labelout(s) = labelin(r);
  labelin(r) = move;
  criterion = [T(N+1,1:K), (-T(N+1,1:K) - 1.0)];
  [margcost, s] = max(criterion);
  count = count + 1;
end
ind = find(abs(labelin) <= K);
blabel = labelin(ind)';
b = T(ind,K+1) .* sign(blabel);
[a, ind] = sort(abs(blabel));
b = b(ind);
b = b+bls;





