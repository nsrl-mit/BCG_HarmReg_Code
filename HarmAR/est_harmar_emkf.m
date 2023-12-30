function [amps_hat,alpha_hat,sigma2_e_iters,sigma2_w_iters,harm_hat,ar_hat,res_hat,ErrCov] = ...
    est_harmar_emkf(w_hat, t, y, vary, N, d, P, A_input, s2w, NextInit, em, ploton)%savename
% Obtain AR and Harmonic Coefficients for Harmonic + Colored Noise Model 
% Input Best Guess of Fundamental Frequency Defining Harmonic
% Uses cyclic descent/E-M to fit data for a given frequency input
% Pavitra Krishnaswamy - 05/2012, Updated 12/2014
%************************************************************************

%% Compute Z for A Given Omega
c = 1; R = 2*d + c;                                                         %   # terms in B_hat
Z = zeros(N,R);                                                             %   y = ZB + v
Z(:,1,:) = ones(N, 1);                                                      %   constant term
%Z(:,2,:) = repmat(t,1,1);                                                  %   take out linear term
nFreq = length(w_hat);
nPerW = floor(length(t)/nFreq);
w = reshape(repmat(w_hat,1,nPerW)',nFreq*nPerW,1);
if length(w)<length(t)
    w = [w;w_hat(end)*ones(length(t)-length(w),1)];
end
for r = 1:d
	a = cos(r*w.*t); Z(:,c+1,:) = a;                                        %   even columns of Z
    b = sin(r*w.*t); Z(:,c+2,:) = b;                                        %   odd columns of Z
    c = c+2;  
end
% Initialize Identity Matrices We Will Need
Ir = eye(R); Ip = eye(P); Irp = eye(R+P);  Ix = eye(N);                     %   identity matrices of size R, P, N
if isempty(A_input)
    A = Ir;                                                                 %   identity harmonic state tx - R X R
else
    A = A_input;                                                            %   non-identity harmonic state tx - R X R
end
Winv = Ix;                                                                  %   initialize inverse AR Covariance
Ipp1 = Ip(1:P-1,:);                                                         %   for AR companion matrix

%% Initialize Cyclic Descent/E-M Loops
count = 1;                                                                  %   count iterations of cyclic descent/E-M
tolem = 1e-2;                                                               %   tolerance for convergence of E-M
tol = 1e-2;
maxiters = 250;                                                             %   max. # iterations for cyclic descent/E-M
keep_loop_running = 1;                                                      %   indicator variable that decides whether to keep loop running

%Initialize AR term as white noise
alpha_hat = zeros(P,1);                                                     %   initialize as white noise
sigma2_e = vary;                                                            %   initialize as white noise (full power)
Psi_Tx = [alpha_hat'; Ipp1];                                                %   update AR state tx (companion/Akaike Matrix) - P X P
sigma2_w = s2w;                                                             %   initialize harmonic amplitude variance for E-M iterations

sigma2_w_iters(count) = sigma2_w;                                           %   for debugging or store only
alpha_hat_iters{count} = alpha_hat;                                         %   for debugging or store only
sigma2_e_iters(count) = sigma2_e;                                           %   for debugging or store only

%% Cyclic Descent Loop/E-M Loop if Running E-M
while(keep_loop_running && (count < maxiters))
	%tic; 
    % Inputs 
    sigma2_w_prev = sigma2_w;                                               %   update sigma2w_prev
    alpha_hat_prev = alpha_hat;                                             %   update alpha_hat_prev
    sigma2_e_prev = sigma2_e;                                               %   update sigma2e_prev
    
    % Initialize State for E-Step: Kalman Recursions
    %updates with iterations automatically due to Winv
    amps_init = (Z'*Winv*Z)\(Z'*Winv*y);                                    %   initialize harmonic amplitudes
    V_init = y-Z*amps_init;                                                 %   initialize AR time series
    x_init = [amps_init; V_init(P:-1:1)];                                   %   augmented state initialization
    %initial state error covariance always use the current sigma2_w and sigma2_e regardless of window or iteration
    Wb = sigma2_w*Ir;                                                       %   state noise covariance for harmonic amplitudes 
    We = 1/100*sigma2_e*Ip; We(1,1) = sigma2_e;                             %   state noise covariance for AR term  --  [sigma2_e*Ip];           
  
    % Bare Minimum E-Step: Run Kalman Filter, Smoother and 1-Step Covariance Algorithms
    [harmar_hat, harm_hat, ar_hat, res_hat, amps_hat, ErrCov.Amps, ar_hat2, ErrCov.AR,...
        ErrCov.Amps_cr, ErrCov.AR_cr] = kalm_filt_smooth(y, Z, A, Psi_Tx, R, P, x_init, Wb, We, Irp, ploton);
    ar_hat      = detrend(ar_hat,'constant');
    %harm_hat   = detrend(harm_hat,'constant');

    % Bare Minimum M-Step: Estimate AR Coefficients Alpha_Hat and AR Covariance
    [arcoeff, sigma2_e] = arburg(ar_hat,P);                                 %   estimate alpha and sigma2e with burg
    alpha_hat = -arcoeff(2:end)';                                           %   update vector of AR coefficients
    Psi_Tx = [alpha_hat'; Ipp1];                                            %   update AR state tx (companion/Akaike Matrix) - P X P
    [Winv, ~, ~] = calcinvcov(ar_hat,P,N,Ix,Ix);                            %   need updated Winv to correct GLS x_init
    %update W_inv with a detrended V so arburg works well
    
    % Full E-M - If EM_FLAG 1
    if em                                  
        % E-Step: Compute Relevant Moments
        [tbb, tb1b1, tbb1] = estep_moments(A, R, N, Wb, amps_hat, x_init, ErrCov);

        % M-Step: Update Values of Sigma2w and Sigma2E
        tw0 = sigma2_w*R-x_init(1:R,1)'*x_init(1:R,1);                      %   initial covariance correction
        sigma2_w = mstep_compute(tbb, tb1b1, tbb1, tw0, R, N);
        
        % Differential with Respect to Last Values
        dsigma2_w = abs((sigma2_w- sigma2_w_prev)/sigma2_w_prev);           %   perc. difference between sigma2_w updates
    else
        % Differential with Respect to Last Values
        dsigma2_w = 1;                                                      %   dummy to keep it greater than tol
    end
    dsigma2_e = abs((sigma2_e -sigma2_e_prev)/sigma2_e_prev);                   %perc. difference between sigma2_e updates
    
    % Decide Whether to Continue Loop Running and Update Counters
    if em
        keep_loop_running = dsigma2_w > tolem || dsigma2_e > tol;
    else
        keep_loop_running = dsigma2_e > tol;
    end
    
    count = count+1;
    sigma2_w_iters(count) = sigma2_w;                                       %   for debugging or store only
    alpha_hat_iters{count} = alpha_hat;                                     %   for debugging or store only
    sigma2_e_iters(count) = sigma2_e;                                       %   for debugging or store only
    %toc;
    
    % Plot to Check Progress
    if ploton     
        figure, set(gcf,'color','white');    
        subplot(2,1,1)
        plot(y); hold on; plot(ar_hat,'r');
        legend('meas','AR'); hold off
        subplot(2,1,2)
        plot(y); hold on; plot(harm_hat,'r')
        legend('meas','Harm'); hold off
        pause(0.5); close;
    end
end

end