% Code plotting pdf from piecewise linear 
clear;
close all;

load Ye_GA_3mix_3_WLS_ue_combined
ntau = 15;

iter  = j;
nmixtures = 3;
lambda_m = nan(iter, nmixtures);
mu_m = nan(iter, nmixtures);
sigma_m = nan(iter, nmixtures);

for j_iter = [1:iter]
    j_iter
    lambdapre  =  (recorder_pc_start_repeat(j_iter,[3*ntau+1:3*ntau+2]));
    mupre  =  (recorder_pc_start_repeat(j_iter,[3*ntau+3:3*ntau+4]));
    sigma  =  (recorder_pc_start_repeat(j_iter,[3*ntau+5:3*ntau+7]));
    
   [lambda,mu] = preprocesslambdamu(lambdapre,mupre); 
   [mu_sorted,mix_index] = sort(mu);
   
   lambda = lambda(mix_index);
   mu = mu(mix_index);
   sigma = sigma(mix_index);
   
   lambda_m(j_iter,:) = lambda;
   mu_m(j_iter,:) = mu;
   sigma_m(j_iter,:) = sigma;
   
end





ParaVector = [mean(lambda_m),mean(mu_m),mean(sigma_m)];
ParaVectortruth = [0.5 0.25 0.25 -3 2 4 1 1 1];



X = [-8:0.01:8];
PDFVector = X*0;
nmixture = 3;
for j_mixture = [1:nmixture]
 PDFVector = PDFVector + ParaVector(j_mixture)*normpdf(X,ParaVector(j_mixture+nmixture),ParaVector(nmixture+nmixture+j_mixture));
end

plot(X,PDFVector,'LineWidth',1.3);

hold on;



PDFVectortruth = X*0;
nmixture = 3;
for j_mixture = [1:nmixture]
 PDFVectortruth = PDFVectortruth + ParaVectortruth(j_mixture)*normpdf(X,ParaVectortruth(j_mixture+nmixture),ParaVectortruth(nmixture+nmixture+j_mixture));
end

plot(X,PDFVectortruth,'r','LineWidth',1.3);
legend('PDF from Average Estimated Measurement Error Parameters','Truth')

ylabel('PDF')
xlabel('x')
axis([-8 8 0 0.6])

print('-dpng','-r0','pdf_GA_truth_1');

