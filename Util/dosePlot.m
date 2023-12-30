function dosePlot(tt,ff,sFull,sHR,sAR,fR,fH,fA,i,N,fRange,cLim,doseStringCell)
%%  Plot raw spectrogram
figure(fR);
a = subplot(N,1,i);
imagesc(tt,ff,10*log10(sFull),cLim);
colormap jet
colorbar;
set(a,'YLim',fRange,'YDir','normal')
xlabel('time (min)','FontSize',10)
ylabel('frequency (Hz)','FontSize',10)
title([doseStringCell;'Data, pre-fit'],'FontSize',12)
%%  Plot harmonic regression fit spectrogram
figure(fH);
b = subplot(N,1,i);
imagesc(tt,ff,10*log10(sHR),cLim);
colormap jet
colorbar;
set(b,'YLim',fRange,'YDir','normal')
xlabel('time (min)','FontSize',10)
ylabel('frequency (Hz)','FontSize',10)
title([doseStringCell;'HR Fit'],'FontSize',12)
%%  Plot autoregressive fit spectrogram
figure(fA);
c = subplot(N,1,i);
imagesc(tt,ff,10*log10(sAR),cLim);
colormap jet
colorbar;
set(c,'YLim',fRange,'YDir','normal')
xlabel('time (min)','FontSize',10)
ylabel('frequency (Hz)','FontSize',10)
title([doseStringCell;'AR Fit'],'FontSize',12)
%%  limit to 30 minutes for methods figures
linkaxes([a b c],'x')
set(a,'XLim',[0 30])