function [lambdasorted, musorted, sigmasorted] = sort_dist(para,ntau)


    lambdapre  =  (para(1,[3*ntau+1 : 3*ntau+2]));
    mupre  =  (para(1,[3*ntau+3 : 3*ntau+4]));
    sigma  =  (para(1,[3*ntau+5 : 3*ntau+7]));
    
    % Sort and save the processed parameters
    [lambda, mu] = preprocesslambdamu(lambdapre,mupre);
    [musorted, I] = sort(mu);
    lambdasorted = lambda(I);
    sigmasorted = sigma(I);
    
    lambdasorted(1,:) = lambdasorted;
    musorted(1,:) = musorted;
    sigmasorted(1,:) = sigmasorted;
end
    
    
