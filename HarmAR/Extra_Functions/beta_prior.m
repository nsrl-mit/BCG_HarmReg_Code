function [Beta_Pri_Inv] = beta_prior(y,w,Fs,d,R)
% Compute Source Covariance for Harmonic Amplitudes
% Pavitra Krishnaswamy - 02/2014
%************************************************************************
    
%% Set Spectral Parameters and Compute Residual Spectrum
NW = 1; K = 2*NW-1; NN = length(y);
params = struct('tapers',[NW K], 'pad', -1, 'Fs', Fs, 'fpass', [0 Fs/2],...
	'movingwin',[NN/Fs NN/Fs],'err',0,'parrun',1);
specres = NW*2*Fs/NN;                                                       %   Hz
[py,fy] = mtspectrumc(y,params);                                            %   uV^2/Hz
res_ind = specres/mean(diff(fy));                                           %   Indices
    
%% Frequencies of Harmonic Peaks
for r = 1:d        
    [~,ind(r)] = min(abs(fy - r*w));
    indres{r} = ind(r)-res_ind:1:ind(r)+res_ind;
end
harm_ind = cell2mat(indres);
    
%% Calculate Source Covariance
Beta_Pri_Inv = zeros(R,R);
c = 3;                                                                      %   set only from 2nd harmonic onward, because 1st harm has slow power
for r = 2:d        
    % indl to indh  defines 2Hz band around rth harmonic
    [~,indl] = min(abs(fy - (r*w - 1)));
    [~,indh] = min(abs(fy - (r*w + 1)));
    ind_band_of_interest = indl:1:indh;
    hyes = ismember(ind_band_of_interest,harm_ind);
    % compute power in bands
    totpow_in_band = py(ind_band_of_interest(hyes));
    bkgdpow_in_band = py(ind_band_of_interest(~hyes));
    desired_harmpow_in_band = mean(totpow_in_band) - mean(bkgdpow_in_band);
    % compute harmonic amplitude prior matrix
    if desired_harmpow_in_band > 0
        Beta_Pri_Inv(c+1,c+1) = 1/desired_harmpow_in_band;
        Beta_Pri_Inv(c+2,c+2) = 1/desired_harmpow_in_band;
    end
    c = c+2;                                                                %   update counter
end
end