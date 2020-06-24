

load MLE_MM_Esample lambdavector_true sigmavector_true muvector_true nmixtures_truth sgrid xgrid sv_true nmixtures l_s l_x
sigmavector_true = [0.5,0.5,0.5]
h = 0.01;
nsample = 10^6;
%% 2) Simulate the data
tau_simu = rand(1,nsample);
% beta0_simu = b0(tau_simu);
beta1_simu = b1_1(tau_simu);
beta2_simu = b2_1(tau_simu);

x2r = rand(1,nsample)+1;


y_ntemp = zeros(1,nsample);
j_sample_begin = 1;
for j_mixture = [1 : nmixtures_truth]
    if j_mixture == nmixtures_truth
        j_sample_end = nsample;
    else
        j_sample_end = min( ceil(sum(lambdavector_true(1:j_mixture))*nsample) , nsample);
    end
    y_ntemp([ j_sample_begin:j_sample_end]) = normrnd(muvector_true(j_mixture),sigmavector_true(j_mixture),1,[(j_sample_end - j_sample_begin +1)]);
    j_sample_begin = j_sample_end+1;
end
[temp , y_n_index] = sort(rand(1,nsample));

y_n = y_ntemp(y_n_index);
y_s = beta1_simu + beta2_simu.*x2r;
y = y_n+y_s;




  for j_s =  [1:l_s]
        s = sgrid(j_s);
        phiep_truth(j_s) = phiep(s,sv_true,nmixtures);
  end
    
  
for j_s = [1:l_s]
     j_s
    s = sgrid(j_s);
    for j_x = [1:l_x]
        x = xgrid(j_x);
        phi_match(j_s,j_x) =  phiY(s,x,h,y,x2r)/phiY(s,x,h,y_n,x*ones(1,nsample));
    end
end





for j_s = [1:l_s]
     j_s
    s = sgrid(j_s);
    for j_x = [1:l_x]
        x = xgrid(j_x);
        phi_match_1(j_s,j_x) =  phiY(s,x,h,y_s,x2r);
    end
end


for j_s = [1:l_s]
   
    s = sgrid(j_s);
    for j_x = [1:l_x]
        x = xgrid(j_x);
        phiT(j_s,j_x) =  phiTrue(s,x);
    end
end
abs(phi_match./phiT)
