function calc_spectra_demo(fnameIn,fNameOut,chan)%,P

for lps = 1:4
    %% Load Data and Initialize
    % Concatenate Across Channels
    for cc = 1:length(chan)
	ch = chan(cc)
    if exist(fnameIn,'file')==2
      load(fnameIn,'raw_hat','ar_hat','harm_hat','noise_hat','time_hat',...
          'starind','stopind','numsegs','Fs'); 
      switch lps 
        case 1 %Raw
          anal_hat = raw_hat;
        case 2 %AR
          anal_hat = ar_hat;
        case 3 %Harmonics
          anal_hat = harm_hat;
        case 4 %Noise 
          anal_hat = noise_hat;
      end
        
      % Block by Block
      for blk = 1:numsegs   
        if isempty(anal_hat{blk}) == 0    
          sig2anal{blk}(:,ch) = anal_hat{blk};
        end
      end
    end
    end
    
    %% Spectral Calculation
    % Initialize
    init = cell(1,numsegs); mtparams = init; ttsd = init; ppsd = init; ffsd = init;  specres = init;
    wbh = waitbar(0,'calculating...');
    % Block by Block
    for blk = 1:numsegs
      if isempty(anal_hat{blk}) == 0 
        % Set Parameters
        dataseg = sig2anal{blk};
        NN = size(dataseg,1);
        resdes = 0.40; NW1 = min(1.5,0.5*ceil(3001*resdes/Fs));
        NW = round(NW1*NN/3001)*(NN > 3001) + NW1*(NN <= 3001); 
        K = 2*NW-1; 
        %bw = 2; NW = 4; K = 2*NW-1;                     %bw means freq resn is 2 Hz
        %[~, N] = dyadlength([1:ceil(Fs*NW/(bw/2))]'); NFFT = 2^N;
        %movingwin = [4 0.25]; ploton = 1;
        params = struct('tapers',[NW K], 'pad', -1, 'Fs', Fs, 'fpass', [0 Fs/2],...
                'movingwin',[NN/Fs NN/Fs],'err',0,'parrun',1);
        
        % Read Outputs from Spectral Calculator
        [pest, test, fest] = mymtspecgram(dataseg,params);
        ppsd{blk} = pest;
        ttsd{blk} = starind(blk)/Fs + test;
        ffsd{blk} = fest;
        mtparams{blk} = params;
        specres{blk} = NW*2*Fs/NN;
        waitbar(blk/numsegs,wbh,'calculating...')
      end
    end
    close(wbh);
    % Store Appropriately
    switch lps 
      case 1 %Raw
        ttsd_r = ttsd; ffsd_r = ffsd; ppsd_r = ppsd;  
        time_hat_r = time_hat; mtparams_r = mtparams; specres_r = specres;
      case 2 %AR
        ttsd_a = ttsd; ffsd_a = ffsd; ppsd_a = ppsd;  
        time_hat_a = time_hat; mtparams_a = mtparams; specres_a = specres;
      case 3 %Harmonics
        ttsd_h = ttsd; ffsd_h = ffsd; ppsd_h = ppsd;  
        time_hat_h = time_hat; mtparams_h = mtparams; specres_h = specres;
      case 4 %Noise
        ttsd_n = ttsd; ffsd_n = ffsd; ppsd_n = ppsd;  
        time_hat_n = time_hat; mtparams_n = mtparams; specres_n = specres;
    end
    clear sig2anal raw_hat ar_hat harm_hat noise_hat time_hat ...
          Fs specres test fest pest params ...
          mot_art ttsd ffsd ppsd mtparams;
end
    
%% Save Spectral Estimates - 1 File for All Channels 
save(fNameOut, 'numsegs', ...
'ttsd_r', 'ffsd_r', 'ppsd_r', 'time_hat_r', 'mtparams_r', 'specres_r',...
'ttsd_h', 'ffsd_h', 'ppsd_h', 'time_hat_h', 'mtparams_h', 'specres_h',...
'ttsd_n', 'ffsd_n', 'ppsd_n', 'time_hat_n', 'mtparams_n', 'specres_n',...
'ttsd_a', 'ffsd_a', 'ppsd_a', 'time_hat_a', 'mtparams_a', 'specres_a'); %'_P',num2str(P),
end