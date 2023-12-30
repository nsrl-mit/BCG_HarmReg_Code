function [likcost, s2e, logdG, S_N] = find_omega(w, ps, t, y, vary, N, d, p)
% Obtain fundamental frequencies for harmonic + colored noise model 
% Uses cyclic descent to fit data and compute concentrated likelihood across possible frequencies
% Pavitra Krishnaswamy - 05/2012
%************************************************************************

%% Compute Z for Each Frequency in W
c = 1; R = 2*d + c;                                                         %   # terms in B_hat
Z = zeros(N,R,ps);                                                          %   y = ZB + v
Z(:,1,:) = ones(N, ps);                                                     %   constant term
%Z(:,2,:) = repmat(t,1,ps);                                                 %   take out any linear trend
for r = 1:d
	a = cos(t*r*w'); 
    Z(:,c+1,:) = a;                                                         %   even columns of Z
    a = sin(t*r*w'); 
    Z(:,c+2,:) = a;                                                         %   odd columns of Z
    c = c+2;  
end

%% Initialize and Run Cyclic Descent for Each Omega
count = ones(ps,1);                                                         %   count iterations of cyclic descent
S_N = zeros(ps,1); logdG = S_N; s2e = S_N;                                  %   initialize variables for costcalc
Ix = diag(ones(1,N));                                                       %   initialize identity matrix
likcost = zeros(ps,1);                                                      %   initialize concentrated likelihood cost

for i = 1:ps
% parfor i = 1:ps
    [arcoeff, s2e_hat] = arburg(y,p);                                       %   estimate alpha and sigma2e as if no harmonics
    alpha_hat = -arcoeff(2:end)';                                           %   initialize vector of AR coefficients
    sigma2e = 0.1*s2e_hat; %0.1*vary;                                       %   initialize white noise variance
    sigma2e_prev = 1.5*sigma2e;                                             %   initialize sigma2e_prev
    Winv = Ix;                                                              %   initialize inverse AR Covariance for this omega
    while(abs(sigma2e-sigma2e_prev) > (1e-4)*abs(sigma2e_prev))
    	sigma2e_prev = sigma2e;                                             %   update sigma2e_prev
        ZZ = Z(:,:,i);                                                      %   Z is NXRXps - if ps is last dimension no need squeeze_fast
        
        %tic;
        amps_hat = (ZZ'*Winv*ZZ)\(ZZ'*Winv*y);                              %   compute amps (RX1) -- + BPInv
        amps_hat(1) = amps_hat(1) + mean(y-ZZ*amps_hat);                    %   correct for any remaining offset in V (typically very small)
        sig_hat = ZZ*amps_hat;                                              %   estimated harmonics 
        V = y - sig_hat;                                                    %   AR series - is this ar_hat + res_hat; mean(V) ~10^(-16)
        [Winv, D, sigma2e]= calcinvcov(V,p,N,Ix,Ix);                        %   update W_inv with a detrended V so arburg works well
        [arcoeff, s2e_hat] = arburg(V,p);                                   %   estimate alpha and sigma2e with burg
        alpha_hat = -arcoeff(2:end)';                                       %   update vector of AR coefficients
        sigma2e = s2e_hat;                                                  %   update white noise variance
        %toc;
        
        S_N(i) = V'*Winv*V;                                                 %   mean square error - diff. from mse of lscov
        logdG(i) = sum(log(D)) + (N-p-1)*log(sigma2e);                      %   log(det(cov(resid AR(p)))) = log(prod([1/(1./([d1 E q .. q])))
        count(i) = count(i)+1;
    end
    s2e(i) = sigma2e;
    likcost(i) = N*log(s2e(i)) + logdG(i) + S_N(i)/s2e(i);                  %   optimal cost(w_hat) = Compute Concentrated logL(w|y)
end

end