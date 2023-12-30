function plot_gof(Estim,GOF,N,R,P,ds,Fs)

%% Variables to Study
resd = Estim.Resid{R,P};                                        %Residual Time Series
normcump = GOF.NCP(R,P);                                        %Residual NCP

%% Normalized Cumulative Periodogram
figure(7), set(gcf,'color','white'); hold on;         
plot(normcump.w, normcump.coef, 'g','LineWidth',4);             %Residual NCP
plot(normcump.w, normcump.theo, 'k');                           %Theoretical NCP    
plot(normcump.w, normcump.theo + normcump.lim,'k--');           %Theoretical NCP Upper Bound
plot(normcump.w, normcump.theo - normcump.lim,'k--');           %Theoretical NCP Lower Bound
[AR.w, AR.cumper, ~, ~] = ncp_compute(Estim.AR{R,P}, Fs, ds);   %AR Estimate NCP
plot(AR.w, AR.cumper,'r','LineWidth',4);                        %Plot AR NCP
[Harm.w, Harm.cumper, ~, ~] = ncp_compute(Estim.Harm{R,P}, Fs, ds); %Harmonic Estimate NCP
plot(Harm.w, Harm.cumper,'b','LineWidth',4);                    %Plot Harmonic NCP
[Meas.w, Meas.cumper, ~, ~] = ncp_compute(Estim.Data{R,P}, Fs, ds); %Measurement NCP
plot(Meas.w, Meas.cumper,'k','LineWidth',4);                    %Plot Measurement NCP
set(gca,'FontSize',14);
xlabel('Frequency (Hz)','FontSize',14); 
ylabel('Normalized Cumulative Periodogram','FontSize',14);
legend('test series','ideal white noise','conf. bounds', 'conf. bounds',...
       'AR','Harmonic','Observation','Location','SouthEast'); 
legend boxoff;
    
%% Q-Q Plot - Compare to Standard Normal
figure(8), set(gcf,'color','white'); qqplot(resd);              
set(gca,'FontSize',14);
title('Q-Q Plot vs. Standard Normal','FontSize',14); 

%% Residual Autocorrelation
%figure(9), set(gcf,'color','white'); zautocorr(resd, 20, 0, 2);%nACFtaps, 0, 2); 
%essentially uses ACF.coef, ACF.lag, ACF.bnd
%figure(10), set(gcf,'color','white'); zparcorr(resd, nACFtaps, 0, 2); 

%% Lagged Scatter Plot
%figure(11), set(gcf,'color','white');
%CS = 5; scatter(resid,circshift(resid,CS);
end