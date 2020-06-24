function result = mn1(r,ParameterDist)
%mn1 defined as function of mixing probability
global nmixtures

% Preprocessed lambda. Lambda=weights of components, only need to specify
% the first nmixtures-1 weights.
lambdapre=(ParameterDist([1:(nmixtures-1)]))';
% Preprocessed mu. mu= mean of components, only need to specify the first nmixtures-1 means.
mupre=(ParameterDist([(nmixtures):(2*nmixtures-2)]))';
%sigma= st.d of each component, need to specify for all nmixtures components.
sigma=ParameterDist([(2*nmixtures-1):(3*nmixtures-2)]);

	z=0;
	%calculate the weight of the nmixtures^th component

    q=1-sum(lambdapre);

	lambda=[lambdapre;q];

	

	if ( q<0 || q==0 ||q>1)
		z=-100000;
		result = z;
        return;
    end
	%mu is the vector of mean of all components.
	mu=[mupre;(0-(mupre'*(lambdapre))/lambda(nmixtures))];

	%calculate density
	for  i = 1:nmixtures
		%if(sigma[i]<0 |sigma[i]==0){return(-100000)}
		z=z+lambda(i)*normpdf((r-mu(i))/sigma(i))/(sigma(i));
    end
	result = z;
