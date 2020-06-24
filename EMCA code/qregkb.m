function [b,se,hse] = qregkb(y,X,q)

b = quantlsf(X,y,q);

n = length(y);
e = y - X*b;
c = (n^(-1/5))*std(e);

t = (15/(16*c))*((1-((e/c).^2)).^2).*(abs((e/c))<1);
J = X'*spdiags(t,0,n,n)*X;
S = X'*spdiags((q-(e<0)).^2,0,n,n)*X;
V = inv(J)*S*inv(J);
se = sqrt(diag(V));

f0 = mean(t);
Vs = (q*(1-q)/(f0^2))*inv(X'*X);
hse = sqrt(diag(Vs));



