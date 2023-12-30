function [FrqOpt, Params, Estim, ErrCov, NextInit] = fit_harmar_emkf(data_ds,data_fs,bounds,...
          R,P,A_input,nW,s2w,omega_input,what_prev,thresh,psize,NextInit,em_flag, ploton)
% Fit Harmonic + Colored (AR) Noise Model to Input Data
% Uses Cyclic Descent/E-M Combination to Fit Data
% Obtain Estimates of Harmonic Frequency + Amplitude, and AR Coeffiences
% Pavitra Krishnaswamy - 05/2012, Updated 12/2014
%************************************************************************
 
%% Read Data
t_fs = data_fs(:,1);                                                        %   full sampling time series
y_fs = data_fs(:,2);                                                        %   full sampling observation series to fit
vary_fs = var(y_fs);                                                        %   variance of observations for full sample
N_fs = size(data_fs, 1);                                                    %   size of full length series
t_ds = data_ds(:,1);                                                        %   downsampling time series
y_ds = data_ds(:,2);                                                        %   downsampling observation series to fit
N_ds = size(data_ds,1);                                                     %   size of downsampled series
wL = floor(N_ds/nW);                                                        %   length of each subwindow (should be n_ds for nW = 1)

%% Set Frequency Search Space
lbound = 2*pi*0.9*bounds(1);
ubound = 2*pi*1.1*bounds(2);
u = linspace(lbound,ubound,psize)';                                         %   span the entire search space uniformly

%% Optimize Cost for Frequency
if omega_input == 0
    w_hat = zeros(nW,1);
    for w = 1:nW
        wRange = (1:wL)+(w-1)*wL;
        vary_ds = var(y_ds(wRange));                                        %   variance of observations for downsampled version
        [cost, s2e, logdG, Sn] = find_omega(u, psize, t_ds(wRange), y_ds(wRange), vary_ds, wL, R, P);                         
        [~, best_freq] = min(cost);                                         %   opt freq is that which minimizes s2e
        x = u(best_freq);                                                   %   optimal best frequency
        w_hat(w) = x/(2*pi);                                                %   convert to Hz from radian frequency
    end
    FrqOpt.lbound = lbound; FrqOpt.ubound = ubound; FrqOpt.u = u; FrqOpt.cost = cost; 
    FrqOpt.s2e = s2e; FrqOpt.logdG = logdG; FrqOpt.Sn = Sn;
else
    w_hat = omega_input; x = w_hat*2*pi;                                    %   use the omega that is input into this function      
    FrqOpt.lbound = lbound; FrqOpt.ubound = ubound; FrqOpt.u = u; FrqOpt.cost = [];
end

%% Check that Optimal Frequency Estimate is in Right Range (not on em trial)
if ~em_flag
    mm = mean(what_prev); 
    if isnan(mm) == 0
        perc_var = (w_hat-mm)/mm*100; 
        if abs(perc_var) > thresh
            w_hat = mm;                                                     %   reset fundamental frequency - its off
            FrqOpt.Misconv = 1;                                             %   flag block for future reference
        else
            FrqOpt.Misconv = 0;
        end
    end
end

%% Optimal Parameter Estimates
[amps,arparams,sigma2_e_iters,sigma2_w_iters,sig_est,ar_est,noise_est,ErrCov] = ...
est_harmar_emkf(2*pi*w_hat, t_fs, y_fs, vary_fs, N_fs, R, P, A_input, s2w, NextInit, em_flag, ploton); 
u_hat = amps(1,:)'; A_hat = amps(2:2:2*R,:)'; B_hat = amps(3:2:(2*R+1),:)'; 
alpha_hat = arparams;

%% Update Initialization for Next Window
NextInit.State = [amps(:,end)' ar_est(end:-1:end-P+1)']';
NextInit.ARParam = alpha_hat;
NextInit.s2e    = sigma2_e_iters(end);
NextInit.propar = sqrt(var(y_fs - sig_est)/vary_fs);
%NextInit.props2e = sqrt(sigma2_e_iters(end)/vary_fs); %residual proportion - used for prev. testing
%NextInit.ErrCov(1:2*R+1,1:2*R+1) = ErrCov.Amps(:,:,end);
%NextInit.ErrCov(2*R+2:2*R+P+1,2*R+2:2*R+P+1) = ErrCov.AR(:,:,end);
%NextInit.ErrCrCov(1:2*R+1,1:2*R+1) = ErrCov.Amps_cr(:,:,end);
%NextInit.ErrCrCov(2*R+2:2*R+P+1,2*R+2:2*R+P+1) = ErrCov.AR_cr(:,:,end);

%% Parameter, Signal and Noise Estimates
Params.Offset = u_hat;                                                      %   Offset
Params.A = A_hat;                                                           %   Cosine Amplitudes
Params.B = B_hat;                                                           %   Sine Amplitues
Params.Freq = w_hat;                                                        %   Fundamental Frequency
Params.Alpha = alpha_hat;                                                   %   AR Coefficients
Params.s2e = sigma2_e_iters(end);                                           %   Residual Noise Variance
Params.s2e_iters = sigma2_e_iters;                                          %   The iterations of s2w values that we went though
Params.s2w = sigma2_w_iters(end);                                           %   The final s2w value used
Params.s2w_iters = sigma2_w_iters;                                          %   The iterations of s2w values that we went though

Estim.Time = t_fs;                                                          %   Time Vector
Estim.Raw = y_fs;                                                           %   Raw Data Vector
Estim.Harm = sig_est;                                                       %   Harmonics Fit
Estim.ARNoise = y_fs - sig_est;                                             %   y-Harmonics Fit
Estim.AR = ar_est;                                                          %   AR Fit
Estim.Noise = noise_est;                                                    %   Residual Noise
end