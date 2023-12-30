function [Estim, HarmAR, AR, GOF] = ...
modelsel_gof(data_ds, data_fs, sub_details, N, harm_only, confLim,nACFtaps, Rst, Ren, Pst, Pen, thresh, ds, Fs)

if harm_only == 1
for R = Rst:Ren
    %% Fit and Estimate Harmonic Only Model
    [~, par_harm, est_harm, likcost] = fit_harm(data_fs, N, sub_details, R);
    w_hat = par_harm.Freq;          Estim.Freq(R) = w_hat;	
    sigma2e = par_harm.Sigma;       Estim.ResVar(R) = sigma2e;
    resd = est_harm.Resid;          Estim.Resid{R} = resd;
    disp(strcat('current model order: R = ', num2str(R)));
    
    %% Model Selection for Joint Harmonics + AR Model: 
    Harm.AIC(R) = N*log(sigma2e) + 2*(2*R+1);
    Harm.AICc(R) = Harm.AIC(R) + 2*(2*R+1)*(2*R+2)/(N-2*R-1-1);
    Harm.BIC(R) = N*log(sigma2e) + (2*R+1)*log(N);
end
end

for R = Rst:Ren
    for P = Pst:Pen
    %% Fit and Estimate Harmonic and AR Contributions
    [~, par_harmar, est_harmar,likcost] = fit_harmar(data_ds, data_fs, sub_details, R, P, [], thresh);
    w_hat = par_harmar.Freq;        Estim.Freq(R,P) = w_hat;
    sigma2e = par_harmar.Sigma;     Estim.ResVar(R,P) = sigma2e;
    harm = est_harmar.Harm;         Estim.Harm{R,P} = harm;
    arseries = (est_harmar.AR)';    Estim.AR{R,P} = arseries;
    resd = est_harmar.Noise;        Estim.Resid{R,P} = resd;
    meas = harm + arseries + resd;  Estim.Data{R,P} =  meas;
    disp(strcat('current model order: R = ', num2str(R), '  P = ', num2str(P)))
    
    %% Model Selection Piece for AR
    AR.AIC(R,P) = N*log(sigma2e) + 2*P;         %N*log(sigma2e/(N-P)) + 2*P;
    AR.AICc(R,P) = AR.AIC(R,P) + 2*P*(P+1)/(N-P-1); %Correct for overfitting if N < 10*(ords^2) - replace P by P+1
    AR.BIC(R,P) = N*log(sigma2e) + P*log(N);    
    AR.FPE(R,P) = sigma2e* (N+P+1)/(N-P-1);     

    %% Model Selection for Joint Harmonics + AR Model: 
    HarmAR.AIC(R,P) = N*log(sigma2e) + 2*(2*R+1+P);
    HarmAR.AICc(R,P) = HarmAR.AIC(R,P) + 2*(2*R+1+P)*(2*R+1+P+1)/(N-2*R-P-1-1);
    HarmAR.BIC(R,P) = N*log(sigma2e) + (2*R+1+P)*log(N);

    %% Test if Residual White or not
    [acf_resid, lbq_resid, ncp_resid] = test_residual(resd, confLim, nACFtaps, P, ds, Fs);
    GOF.ACF(R,P) = acf_resid;
    GOF.LBQ(R,P) = lbq_resid;
    GOF.NCP(R,P) = ncp_resid;
    Estim.LBQ_Res(R,P) = lbq_resid.H;
    end
end

end