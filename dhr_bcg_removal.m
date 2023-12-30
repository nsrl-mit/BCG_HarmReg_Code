function [RawEEG, CleanEEG, ArtBCG] = dhr_bcg_removal(eegdata, chan, fs, heartbeat, ploton, dataSavePath, saveon, plot_duration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%outputs clean_s and artifact_same size as eeg data, same sampling frequency and in continuous formats 
%cleans iteratively channel by channel
%Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu), February 1, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Set Log File and Saving DetailsPaths for Data, Dynamic Harmonic Regression, and Chronux and Plotting Code Streams
saveFiletag             = 'Results_Chan';                                   %   Channel wise save filename
logFile                 = [dataSavePath, 'fitLog.txt'];                     %   Log file to save computing times and diagnostics
efh = fopen(logFile,'a');
fprintf(efh,'%s\t\t\t\t%s\t\t%s\t\t%s\n','Date and Time', 'Channel #', 'Estimation Status','Run Time (Minutes)');
fclose(efh);
        
%%  Set Global Params
% Data Processing Parameters: Channels and Time Windows
dsFactor                = 2;                                                %   downsampling factor, used in frequency estimation
tStart                  = 0;                                                %   starting time  (s)
tTotal                  = 120;                                              %   total time (s)
tWindow                 = 3;                                                %   duration of moving window (s)
wLength                 = fs*tWindow;                                       %   number of samples per window
nWindow                 = tTotal/tWindow;                                   %   number of windows considered
% Dynamic Harmonic Regression Model and Estimation Parameters
R                       = 18;                                               %   order of harmonic regression model component
P                       = 3;                                                %   order of autoregressive model component
A                       = eye(2*R+1);                                       %   harmonic state evolution transform
omega                   = 0;                                                %   initial frequency estimate (if 0, the model estimates at each window)                
if exist('Events_Heart_Beats.xml','file') && heartbeat == 1
	bounds              = [];                                               %   read QRS peaks from Events_Heart_Beats.xml
else
    bounds              = [50 100]/60;                                      %   bounds of frequency estimation: naive set to 50-100 BPM
end
pSize                   = 100;                                              %   granularity of frequency estimates (range(bounds)/pSize)
thresh                  = 10;                                               %   thresholded change of fundamental frequency estimates between windows
% EM Parameters
emPeriod                = 60;                                               %   time between EM implementations (s) (set to 0 if EM is disabled)
emWindow                = 15;                                               %   window length used for EM implementation (s) - keep it a multiple of tWindow
s2wFrac                 = 1e-7;                                             %   scaling factor to determine initial s2w
pwrPeriod               = 60;                                               %   amount of time used to get a power baseline for s2w estimation (seconds)
% Initialization
RawEEG                  = zeros(size(eegdata'));                            %   initialize output series
CleanEEG                = zeros(size(eegdata'));                            %   initialize output series
ArtBCG                  = zeros(size(eegdata'));                            %   initialize output serie

%%  Run Dynamic Harmonic Regression for Each Individual Channel
for k = 1:size(eegdata, 2);                                                 %   All Channels 
    time_series         = [0:1:length(eegdata(:,k))-1]/fs;                  %   Time Series
    EEG_data            = [time_series; eegdata(:,k)']';                    %   EEG Series
   
    if ismember(k, chan)
      try 
        %   Run Estimation
        % Initialize BCG Power/Signal Power Ratio
        s2w(k)          = s2wFrac*sum(EEG_data(1:round(pwrPeriod*fs),2).^2)/(pwrPeriod);
    
        % Call Estimation for this Channel 
        ctic            = tic;
        [EEG_Series, FrqOpt, Params, Estim, ErrCov, NextInit] = fitWinEMKF(EEG_data, floor(fs), bounds, pSize, thresh, dsFactor,...
            tWindow, R, P, A, s2w(k), omega, emPeriod, emWindow, ploton, plot_duration); 
        compTime        = toc(ctic);
        
        % Save Outputs
        RawEEG(k, :)    = EEG_Series.Raw';                                  %   Raw EEG Series
        CleanEEG(k, :)  = EEG_Series.Clean';                                %   Clean EEG Series
        ArtBCG(k, :)    = EEG_Series.BCG';                                  %   BCG Artifact Series       
        
        % Save Detailed Results for this Channel
        if saveon
        disp(strcat('now saving results for channel...', num2str(k), '...'));
        saveFileName = strcat(saveFiletag, num2str(k),'.mat');              %   Channel Wise Saving
        save([dataSavePath saveFileName],'FrqOpt','Params','Estim','ErrCov','NextInit','fs','-v7.3');
        end
        
        % Write Successful Fit to Log File
        efh = fopen(logFile,'a');
        fprintf(efh,'%s\t\t%s\t\t%s\t\t\t%f\n',char(datetime()),strcat('channel...',num2str(k)),'est. successful',floor(compTime/60));
      catch fittingError                                                     %   Log if any Errors in Estimation
        % Write fitting error to log file
        efh = fopen(logFile,'a');
        fprintf(efh,'%s\t\t%s\t\t%s\t\t\t%f\n',char(datetime()),strcat('channel...',num2str(k)),fittingError.message, floor(toc(ctic)/60));
        fclose(efh);
      end
    else                                                                    %   Return the Data as is
        RawEEG(k, :)    = EEG_data(:,2);                                    %   Raw EEG Series
        CleanEEG(k, :)  = EEG_data(:,2);                                    %   Clean EEG Series
        ArtBCG(k, :)    = zeros(size(EEG_data(:,2)));                       %   BCG Artifact Series
    end
end

end