% Code plotting pdf from piecewise linear 
clear;
close all;

load Ye_GA_3mix_3_WLS_ue_combined

iter = 500;
nmixture = 3;

X = [-8:0.01:8];

PDF_matrix = zeros(iter,length(X));
pl_dist_unique = recorder_pc_start_repeat(:,[46:52]);
for j_iter = [1:iter]
    lambdapre  =  pl_dist_unique(j_iter,[1,2]);
    mupre  =  pl_dist_unique(j_iter,[3:4]);
    sigma  =  pl_dist_unique(j_iter,[5:7]);
    
   [lambda,mu] = preprocesslambdamu(lambdapre,mupre); 
   [mu_sorted,mix_index] = sort(mu);
   
   lambda = lambda(mix_index);
   mu = mu(mix_index);
   sigma = sigma(mix_index);
   

for j_mixture = [1:nmixture]
 PDF_matrix(j_iter,:) =  PDF_matrix(j_iter,:) + lambda(j_mixture)*normpdf(X,mu(j_mixture),sigma(j_mixture));
end
   

end

ParaVectortruth = [0.5 0.25 0.25 -3 2 4 1 1 1];




plot(X,mean(PDF_matrix),'LineWidth',1.3);

hold on;



PDFVectortruth = X*0;
nmixture = 3;
for j_mixture = [1:nmixture]
 PDFVectortruth = PDFVectortruth + ParaVectortruth(j_mixture)*normpdf(X,ParaVectortruth(j_mixture+nmixture),ParaVectortruth(nmixture+nmixture+j_mixture));
end

plot(X,PDFVectortruth,'r','LineWidth',1.3);
legend('Average PDF','Truth')

v1 = (quantile(PDF_matrix,0.025));
v2 = (quantile(PDF_matrix,0.975));
plot(X,v1,'b--','LineWidth',1);
plot(X,v2,'b--','LineWidth',1);


output_matrix = [X; PDFVectortruth;mean(PDF_matrix);  v1; v2]';

csvwrite('pdf.csv', output_matrix);

ylabel('PDF')
xlabel('x')
axis([-8 8 0 0.3])

print('-dpng','-r0','pdf_GA_truth_1');


