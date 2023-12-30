function plot_modsel(Estim, HarmAR, AR, Rst, Ren, Pst, Pen)
% Plot model selection results for harmonic regression in colored noise
% Pavitra Krishnaswamy - 08/2013
%************************************************************************

%% Model Selection
%2D Joint BIC Plot
figure(1), set(gcf,'color','white'); imagesc(Pst:Pen,Rst:Ren,HarmAR.BIC(Rst:Ren,Pst:Pen));       
axis xy; set(gca,'FontSize',14); axis([Pst Pen Rst Ren]); colorbar; 
xlabel('AR Order','FontSize',14); ylabel('Harmonic Order','FontSize',14); 
title('Joint BIC vs. Model Orders','FontSize',14); 

%figure, set(gcf,'color','white'); plot(Harm.BIC); 
%xlabel('Harmonic Order'); title('Purely Harmonic Model')  
%descending graph -- as AR has lots of power

%Joint BIC vs. Harm Order
figure(2), set(gcf,'color','white'); plot(HarmAR.BIC);          
[~,b] = min(HarmAR.BIC(Rst:Ren,Pst:Pen),[], 1); h = Rst+round(mean(b))-1 %Find Best Harmonic Order
xlim([Rst Ren]); set(gca,'FontSize',14); 
xlabel('Harmonic Order','FontSize',14); ylabel('Joint Bayesian Information Criterion','FontSize',14);
title('Harmonic + AR Model','FontSize',14);

%Joint BIC vs. AR Order for Best Harmonic Model
figure(3), set(gcf,'color','white'); plot(HarmAR.BIC(h,:));     
xlim([Pst Pen]); set(gca,'FontSize',14); 
xlabel('AR Order','FontSize',14); ylabel('Joint Bayesian Information Criterion','FontSize',14); 
title('With Best Harmonic Model, Choose AR Model','FontSize',14);

%AR Only BIC for Best Harmonic Model
figure(4), set(gcf,'color','white'); plot(AR.BIC(h,:));         
xlim([Pst Pen]); set(gca,'FontSize',14);
xlabel('AR Order','FontSize',14); ylabel('AR Only Bayesian Information Criterion','FontSize',14)
title('With Best Harmonic Model, Choose AR Model','FontSize',14);

%Influence of Model Order on LBQ Test
figure(5), set(gcf,'color','white'); imagesc(Estim.LBQ_Res);    
axis xy; axis([Pst Pen Rst Ren]); set(gca,'FontSize',14); colorbar; 
xlabel('AR Order','FontSize',14); ylabel('Harmonic Order','FontSize',14); 
title('GOF Across Models','FontSize',14); 
%null is no residual autocorrelation - plot = 0 means failure to reject null

%Influence of Model Order on Frequency
figure(6), set(gcf,'color','white');  imagesc(Estim.Freq);      
axis xy; axis([Pst Pen Rst Ren]); set(gca,'FontSize',14); colorbar; 
xlabel('AR Order','FontSize',14); ylabel('Harmonic Order','FontSize',14); 
title('Opt Freq. Vs. Model','FontSize',14);

end