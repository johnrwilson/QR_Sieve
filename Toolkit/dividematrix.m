function [ begin_index,end_index,index_length ] = dividematrix( nsample,nslice )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dividematrix
% Divide one marix into multiple slices
%
% Errors in the Dependent Variable of Quantile Regression Models
%
% Jerry Hausman, Haoyang Liu, Ye Luo, Christopher Palmer 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


neachsample = floor(nsample/nslice);


begin_index = 1+[0:(nslice-1)]*neachsample;
end_index = [1:(nslice)]*neachsample;
end_index(end) = nsample;
index_length = end_index - begin_index+1;
end


