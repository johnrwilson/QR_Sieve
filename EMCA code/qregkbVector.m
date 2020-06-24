function [fit,sev,hsev] = qregkbVector(y,X,lp)
lv = length(lp);
[temp,lx] = size(X);
fit = NaN(lv,lx);
sev = fit;
hsev = fit;
for jv = [1:lv]
    
    [q  ] = lp(jv);
    [b,se,hse] = qregkb(y,X,q);
fit(jv,:) = b;
sev(jv,:) = se;
hsev(jv,:) = hse;
    
end
