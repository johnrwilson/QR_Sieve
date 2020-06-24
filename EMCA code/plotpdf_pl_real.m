% Code plotting pdf from piecewise linear real data result
% 9/22/2018
clear;
close all;


% Part A) Load data
% 80 White 15 tau
 load RData_pl_Angrist80_sub1_ntau15_white1_ue1
% 90 White 15 tau
%load RData_pl_Angrist90_sub1_ntau15_white1_ue1
% 00 White 15 tau
% load RData_pl_Angrist00_sub1_ntau15_white1_ue1
% 10 White 15 tau
% load RData_pl_Jacob10_sub1_ntau15_white1_ue1


% Part B) The PDF result
X = [-8:0.01:8];
PDFVector = X*0;

lambdapre  =  (recorder_pc_start_repeat(1,[4*ntau+1:4*ntau+2]));
mupre  =  (recorder_pc_start_repeat(1,[4*ntau+3:4*ntau+4]));
sigma  =  (recorder_pc_start_repeat(1,[4*ntau+5:4*ntau+7]));

[lambda,mu] = preprocesslambdamu(lambdapre,mupre);

ParaVector = [lambda, mu, sigma];
return;
ParaVector
for j_mixture = [1:3]
    PDFVector = PDFVector + ParaVector(j_mixture)*normpdf(X,ParaVector(j_mixture+3),ParaVector(6+j_mixture));
end

% Part C) Plotting

plot(X,PDFVector,'LineWidth',1.3);
%legend('PDF from Estimated Measurement Error Parameters','Truth')
ylabel('PDF')
xlabel('x')
title((sprintf('Data: %s',postfix)))

print('-dpng','-r0',sprintf('PDF_%s',postfix));


