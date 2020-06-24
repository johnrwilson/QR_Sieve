function result = mnnoprevector(r_vector,ParameterDist,nmixtures)
%mn1 defined as function of mixing probability

% Preprocessed lambda. Lambda=weights of components, only need to specify
% the first nmixtures-1 weights.
lambda=(ParameterDist([1:(nmixtures)]))';
% Preprocessed mu. mu= mean of components, only need to specify the first nmixtures-1 means.
mu=(ParameterDist([(nmixtures+1):(2*nmixtures)]))';
%sigma= st.d of each component, need to specify for all nmixtures components.
sigma=ParameterDist([(2*nmixtures+1):(3*nmixtures)]);


z = zeros(1,length(r_vector));
	%calculate density
	for  i = 1:nmixtures
		%if(sigma[i]<0 |sigma[i]==0){return(-100000)}
		z=z+lambda(i)*normpdf((r_vector-mu(i))/sigma(i))/(sigma(i));
    end
	result = z;
