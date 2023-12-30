function [b,E] = arburg_parvec(x, p)   
%   ARBURG2 Vectorized AR parameter estimation via Burg method
%   x is NXps, p is a number (model order of the AR system)
%   a{m} is polynomial A corresponding to the AR parametric signal model 
%   estimate of vector X using Burg's method, b{m} has fn(coefficients)!
%   E is final prediction error E (estimated white noise variance)
%   b{m} has size psXm and m can range from 1 to p, E is mXps

%   See also PBURG, ARMCOV, ARCOV, ARYULE, LPC, PRONY.
%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
%              Prentice-Hall, 1988, Chapter 7
%        S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
%              Macmillan, 1988, Chapter 5
%   Author(s): D. Orofino and R. Losada, Modified by Pavitra Krishnaswamy
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.12.4.6 $  $Date: 2010/12/06 00:01:12; 2012/06/05 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check the input data type. Single precision is not supported.
narginchk(2,2)
try
    chkinputdatatype(x,p);
catch ME
    throwAsCaller(ME);
end

validateattributes(x,{'numeric'},{'nonempty','finite', 'vector'},'arburg','X');
validateattributes(p,{'numeric'},{'positive','integer','scalar'},'arburg','ORDER');
if issparse(x),
   error(message('signal:arburg:Sparse'))
end
if numel(x) < p+1
    error(message('signal:arburg:InvalidDimension', p + 1));
end

%% Initialize
x  = x(:);                                                                  %   make x a column vector
N  = length(x);                                                             %   #time points in each residual vector
ef = x;                                                                     %   Initial forward term - size NXps
eb = x;                                                                     %   Initial back term - size NXps
b = cell(1,p);                                                              %   Preallocate for speed
a = 1;                                                                      %   Initial parameter - size 1X1
E = diag(x'*x)./N;                                                          %   Initial error - size 1X1
k = zeros(p,1);                                                             %   Preallocate 'k' for speed.

%% Burg Algorithm Recursions
for m=1:p
   % Calculate the next order reflection (parcor) coefficient
   efp = ef(2:end);                                                         %   n-1X1
   ebp = eb(1:end-1);                                                       %   n-1X1
   num = -2.*ebp'*efp;                                                      %   1X1
   den = efp'*efp+ebp'*ebp;                                                 %   1X1
   kk = num./den; k(m) = kk;                                                %   1X1

   % Update the forward and backward prediction errors
   ef = efp + kk*ebp;                                                       %   n-1X1
   eb = ebp + kk*efp;                                                       %   n-1X1
   
   % Update the AR coeff.
   a=[a;0] + kk*[0;conj(flipud(a))];                                        %   length m + 1 - grows by 1 each iteration
   coef = -a(2:end)'; coef = -fliplr(coef); b{m} = coef;                    %   Z to t-domain & flips
   
   % Update the prediction error
   E(m+1) = (1 - kk.*kk)*E(m);                                              %   1X1
end

%% Assign Outputs
% a{m}(n,r) for nth population of inputs, AR(r)
% b{m}(n,r) for nth populatin of inputs, -fliplr(AR(r-1))
E = E(2:end);                                                               %   r^th element for error of fitting to AR(r)

end