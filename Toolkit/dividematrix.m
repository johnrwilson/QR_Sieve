function [ begin_index,end_index,index_length ] = dividematrix( nsample,nslice )

neachsample = floor(nsample/nslice);


begin_index = 1+[0:(nslice-1)]*neachsample;
end_index = [1:(nslice)]*neachsample;
end_index(end) = nsample;
index_length = end_index - begin_index+1;
end


