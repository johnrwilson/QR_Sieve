function [mu] = preprocesslambdamu1(lambda,mupre)


lambdapre = lambda([1:(end-1)]);
lambda3 = lambda(end);
% Process mu
if length(mupre)>=1
mu=[mupre,(0-(mupre*(lambdapre'))/(lambda3))];
return;
end
if isempty(mupre)
    mu = 0;
end
    

