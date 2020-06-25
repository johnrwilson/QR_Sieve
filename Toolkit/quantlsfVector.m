function [fit] = quantlsfVector(X,y,lp)
lv = length(lp);
[temp,lx] = size(X);
fit = NaN(lv,lx);
for jv = [1:lv] 
    [q  ] = lp(jv);
    [b] = quantlsf(X,y,q);
    fit(jv,:) = b;
end
