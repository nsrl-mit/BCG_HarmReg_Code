function [seq_fit, seq_res] = argen(coefs, seq)
%seq  = Vres'; coefs = alpha_hat';
% Generate an AR sequence of length N based on the specified AR
% coefficients and using the first P samples of given sequence as
% seeds for forward prediction, then use backward prediction to calculate
% the first P samples of the predicted sequence
% seq_fit = ARgen(coefs, seq)
% Inputs
%   coefs       AR coefficients (time-domain lag form). (Px1)
%   seq         The sequence to which the AR model has been fit. (1xN)
% Outputs
%   seq_fit     The predicted sequence. (1xN)
%   seq_res     The residuals, i.e. the sample-wise prediction error. (1xN) 
%               seq_res = seq - seq_fit
% Wasim Malik - 12/18/2008
% Modified by Pavitra Krishnaswamy - 05/2012
%--------------------------------------------------------------------------

%% Set variables
N = length(seq);
P = length(coefs);
    
%% Prepare 
fcoefs = flipud(coefs);                                                     %   flip the sequence of coefs
%size(fcoefs)
ind = 1 : (N-P);                                                            %   predictors (previous P values)
mat = [];

%% Forward prediction of elements P+1:end
for i = 0 : P-1                                                             %   make predictor matrix
	mat = [mat; seq(ind+i)];
end
seq_fit = (mat.' * fcoefs).';                                               %   do the fit estimation

%% Back prediction of elements 1:P
for p = 1 : P
	q = P - p + 1;                                                          %   reverse index
    a = seq(q+1 : q+P);                                                     %   read in estimates for this segment
    b = coefs;                                                              %   estimated coefficients
    bk(q) = a * b;                                                          %   back predicted elements
end

%% Put it all together
seq_fit = [bk seq_fit];
seq_res = seq - seq_fit;                                                    %   the noise (error) from the model fit
end