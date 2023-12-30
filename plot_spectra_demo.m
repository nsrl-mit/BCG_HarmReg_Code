function plot_spectra_demo(fIn, chan, dbs, cax)%,P)
  %% Initialize
  st_time_r = 0; st_time_a = 0; st_time_h = 0; st_time_n = 0; st_time_b = 0;
  smp = 3;          %spectral resolution divided by frequency bin size
  if dbs == 0       %better to use dB when need to reduce dynamic range (but we dont need that)
      colorunits = 'PSD (uV^2/Hz)';
  else
      colorunits = 'PSD (dB)';
  end

  %% Load Files
  fname = fIn; %'_P',num2str(P),
  if exist(fname,'file')==2
  load(fname,'numsegs',...
  'time_hat_r','ffsd_r','ppsd_r','time_hat_a','ffsd_a','ppsd_a',...
  'time_hat_h','ffsd_h','ppsd_h','time_hat_n','ffsd_n','ppsd_n');
  time_hat_b = time_hat_n;
  for cc = 1:length(chan)
      %% Setup Channel Specific Parameters
      ch = chan(cc); disp(strcat('channel: ', num2str(ch)));
      sname{cc} = {strcat('eegsp_raw_', 'ch',num2str(ch)),...
                   strcat('eegsp_ar_', 'ch',num2str(ch)),...
                   strcat('eegsp_harm_', 'ch',num2str(ch)),...
                   strcat('eegsp_resid_', 'ch',num2str(ch))};

      %% Plot Raw Spectra
      figure(ch), set(gcf,'color','white'); hold on;
      keyboard
      for blk = 1:numsegs
      if isempty(ppsd_r{blk}) == 0
        if dbs == 1
          imagesc((st_time_r + time_hat_r{blk})/60, ffsd_r{blk}', 10*log10(ppsd_r{blk}(ch,:)));
        else
          imagesc((st_time_r + time_hat_r{blk})/60, ffsd_r{blk}', ppsd_r{blk}(ch,:));
        end
        hold on; axis xy;
        ttr = time_hat_r{blk};
      end
      end
      figure(ch), axis([st_time_r (st_time_r + ttr(end))/60 0.5 25]); caxis(cax); colorbar;
      title(strcat(' Channel', num2str(ch),' Raw Data'),'FontSize',14); hold off;
      xlabel('Time (Minutes)','FontSize',14);          % x-axis label
      ylabel('Frequency (Hz)','FontSize',14);          % y-axis label
      ylabel(colorbar, colorunits,'FontSize',14);      % colorbar label
      set(gca,'FontSize',14);                          % change the font for the ticks

      %% Plot AR Fit Spectra
      % smoothing parameter is derived as specres_a{blk}./mean(diff(ffsd_a{blk}))
      figure(3*ch), set(gcf,'color','white'); hold on;
      for blk = 1:numsegs
      if isempty(ppsd_a{blk}) == 0
        if dbs == 1
          imagesc((st_time_a + time_hat_a{blk})/60, ffsd_a{blk}', 10*log10(ppsd_a{blk}(ch,:)));
          %imagesc((st_time_a + time_hat_a{blk})/60, ffsd_a{blk}', 10*log10(ppsd_a{blk}(:,ch)));
        else
          imagesc((st_time_a + time_hat_a{blk})/60, ffsd_a{blk}', smooth(ppsd_a{blk}(ch,:),smp));
          %imagesc((st_time_a + time_hat_a{blk})/60, ffsd_a{blk}', ppsd_a{blk}(:,ch));
        end
        hold on; axis xy;
        tta = time_hat_a{blk};
      end
      end
      figure(3*ch), axis([st_time_a (st_time_a + tta(end))/60 0.5 25]); caxis(cax); colorbar;
      title(strcat(' Channel', num2str(ch),' AR Fit'),'FontSize',14); hold off;
      xlabel('Time (Minutes)','FontSize',14);          % x-axis label
      ylabel('Frequency (Hz)','FontSize',14);          % y-axis label
      ylabel(colorbar, colorunits,'FontSize',14);      % colorbar label
      set(gca,'FontSize',14);                          % change the font for the ticks

      %% Plot Harmonic Fit Spectra
      figure(5*ch), set(gcf,'color','white'); hold on;
      for blk = 1:numsegs
      if isempty(ppsd_h{blk}) == 0
        if dbs == 1
          imagesc((st_time_h + time_hat_h{blk})/60, ffsd_h{blk}', 10*log10(ppsd_h{blk}(ch,:)));
        else
          imagesc((st_time_h + time_hat_h{blk})/60, ffsd_h{blk}', ppsd_h{blk}(ch,:));
        end
        hold on; axis xy;
        tth = time_hat_h{blk};
      end
      end
      figure(5*ch), axis([st_time_h (st_time_h + tth(end))/60 0.5 25]); caxis(cax); colorbar;
      title(strcat(' Channel', num2str(ch), ' Harmonics'),'FontSize',14); hold off;
      xlabel('Time (Minutes)','FontSize',14);          % x-axis label
      ylabel('Frequency (Hz)','FontSize',14);          % y-axis label
      ylabel(colorbar, colorunits,'FontSize',14);      % colorbar label
      set(gca,'FontSize',14);                          % change the font for the ticks

      %% Plot Noise Spectra
      figure(7*ch), set(gcf,'color','white'); hold on;
      for blk = 1:numsegs
      if isempty(ppsd_n{blk}) == 0
        if dbs == 1
          imagesc((st_time_n + time_hat_n{blk})/60, ffsd_n{blk}', 10*log10(ppsd_n{blk}(ch,:)));
        else
          imagesc((st_time_n + time_hat_n{blk})/60, ffsd_n{blk}', ppsd_n{blk}(ch,:));
        end
        hold on; axis xy;
        ttn = time_hat_n{blk};
      end
      end
      figure(7*ch), axis([st_time_n (st_time_n + ttn(end))/60 0.5 25]); caxis(cax); colorbar;
      title(strcat(' Channel', num2str(ch), ' Residual'),'FontSize',14); hold off;
      xlabel('Time (Minutes)','FontSize',14);          % x-axis label
      ylabel('Frequency (Hz)','FontSize',14);          % y-axis label
      ylabel(colorbar, colorunits,'FontSize',14);      % colorbar label
      set(gca,'FontSize',14);                          % change the font for the ticks
  end
  else
      disp('file does not exist')
  end

  %% Save Figures
  for cc = 1:length(chan)
  fig_ind = 1:2:7;
  for i = 1:length(fig_ind)
    savename = sname{cc}{i};                           % strcat(,'_P',num2str(P));
    figure(fig_ind(i)*ch), h = gcf;
    %saveas(h, strcat(savename,'.jpg'));               % jpg
    print(h,'-dpdf',savename);                         % pdf
    %print(h,'-deps',savename);                        % -deps for greyscale
    %print(h,'-depsc',savename);                       % -depsc for color
  end
  end
end
