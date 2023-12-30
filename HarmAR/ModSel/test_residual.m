function [ACF, LBQ, NCP] = test_residual(ARres, confLim, nACFtaps, p, ds, Fs)
% Test whether the given sequence has white noise characteristics
% Wasim Malik - 07/2008
% Modified by Pavitra Krishnaswamy - 08/2012
%Outputs described below
%ACF.coef is sample autocorrelation function of t
%ACF.lag is vector of lags corresponding to acf
%ACF.bnd indicates approximate upper and lower confidence bounds
%LBQ.H is vector of boolean decisions for tests - h = 0 fail to reject null
%LBQ.pVal is vector of p-values of the test statistics
%LBQ.Qstat is vector of test statistics 
%LBQ.critVal is vector of critical values for the tests
%NCP.w is frequency axis
%NCP.coef is cumulative periodogram values
%NCP.lim is confidence limits around NCP.coef
%NCP.theo is for theoretical white noise
%************************************************************************

%% Set Bounds based on Confidence Limits
switch confLim          % Set coefficient K_e for the AR redisual test based on normalized cumulative periodogram [pp. 324, 1]
	case 99             % 99% confidence
    	K_e = 1.63;
        nsd = 2.54;    
	case 95             % 95% confidence
    	K_e = 1.36;
        nsd = 1.96;     
	case 90             % 90% confidence
    	K_e = 1.22;
        nsd = 1.63;     
	case 75             % 75% confidence
    	K_e = 1.02;
        nsd = 0.85;     
    otherwise
        disp('ARfit: wrong confidence interval specification');
end

%% Autocorrelation test
len_res = length(ARres);
[ACF.coef, ACF.lag, ACF.bnd] = zautocorr(ARres, nACFtaps, 0, nsd);

%% Modified LBQ test for lack of fit. See eq. (8.2.3) in [pp. 314, 1]
nlags = round(log(len_res)); %max(p+1,round(log(len_res))); %nACFtaps
[LBQ.H, LBQ.pVal, LBQ.Qstat, LBQ.critVal] = zlbqtest(ARres, 'lags', nlags, 'alpha',(1 - confLim/100),'dof',max(nlags-p,1));

%% Periodogram test
[NCP.w, NCP.coef, len_res, len_pxx] = ncp_compute(ARres, Fs, ds);
if mod(len_res, 2)      % odd
	q = (len_res - 1) / 2;
else                    % even
	q = (len_res - 2) / 2;
end
NCP.lim = K_e/sqrt(q); 
NCP.theo = linspace(0, 1, len_pxx); % bounds 

end  