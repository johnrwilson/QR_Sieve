function [lambda,mu,lambda3] = preprocesslambdamu(lambdapre,mupre)

% Process lambda
lambda3 = 1-sum(lambdapre);
lambda = [lambdapre,lambda3];

% Process mu
if length(mupre)>=1
mu=[mupre,(0-(mupre*(lambdapre'))/(lambda3))];
return;
end
if isempty(mupre)
    mu = 0;
end
    

