%% Changes the definition of parameters in gradl_CDF_Lei to avoid constraints. 
% 2/3/2018
function [cdf_lf, cdf_g] = gradl_search_Lei(para_raw, taugrid, nmixtures, Y, X)
    para_Pre = para_raw;
    para_Pre(end-2:end) = exp(para_Pre(end-2:end)-1);
    b1 = para_Pre(end-6);
    b2 = para_Pre(end-5);
    w1 = exp(b1)/(exp(b1) + exp(b2) + 1);
    w2 = exp(b2)/(exp(b1) + exp(b2) + 1);
    para_Pre(end-6) = w1;
    para_Pre(end-5) = w2;
   
    [cdf_lf, cdf_g_cdf] = gradl_CDF_Lei(para_Pre, taugrid, nmixtures, Y, X);
    
    cdf_g = cdf_g_cdf;
    cdf_g(end-2:end) = exp(para_raw(end-2:end)-1).*cdf_g_cdf(end-2:end);
    cdf_g(end-6) = cdf_g_cdf(end-6)*exp(b1)*(exp(b2)+1)/((exp(b1) + exp(b2) + 1)^2) - cdf_g_cdf(end-5)*exp(b1)*exp(b2)/((exp(b1) + exp(b2) + 1)^2);
    cdf_g(end-5) = cdf_g_cdf(end-5)*exp(b2)*(exp(b1)+1)/((exp(b1) + exp(b2) + 1)^2) - cdf_g_cdf(end-6)*exp(b2)*exp(b1)/((exp(b1) + exp(b2) + 1)^2);
    
end
