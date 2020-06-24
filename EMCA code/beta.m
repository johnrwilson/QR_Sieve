
% b(t)
function value = beta(theta,tgrid,t,dtau)
    t = t';
    tgrid = repmat(tgrid, length(t), 1);
% 	tgrid is the grid of taus
% 	tvector is the vector of taus that we want to evaluate
% 	theta is the vector of parameters
% 	dtau is step size
	[tbelow,tbelow_index]=max((tgrid<t).*tgrid, [], 2);
    thetamatrix = repmat(theta, length(t), 1);
    indexmat = ones(length(t), length(theta));
    indexmat = (cumsum(indexmat,2)<=tbelow_index);
    indexmat(:,1) = 0 ;
    constantterm = sum(thetamatrix.*indexmat,2)*dtau;
	value = theta(1)+constantterm+theta(tbelow_index+1)'.*(t-tbelow);
end
