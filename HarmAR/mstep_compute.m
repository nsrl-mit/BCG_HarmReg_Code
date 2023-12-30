function sigma2_w = mstep_compute(tbb, tb1b1, tbb1, tw0, R, N)
% Use Moments from E-Step to Maximize Complete Data Likelihood in Unknown Parameters
% Pavitra Krishnaswamy 12/2014
%************************************************************************

% sigma2_w  = 1/((N+1)*R)*(sum(tbb) + sum(tb1b1) - 2*sum(tbb1) + tw0);
sigma2_w = 1/((N+1)*R)*(sum(tbb) + sum(tb1b1) - 2*sum(tbb1));

end