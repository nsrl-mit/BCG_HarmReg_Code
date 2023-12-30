function [EEG_Series, FrqOpt,Params,Estim,ErrCov,NextInit,N] = fitWinEMKF(data,fs,bounds,pSize,thresh,dsFactor,...
         windowTime,R,P,A_input,sigma2w,omega,emPeriod,emWindow, ploton, plot_duration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data: [time points; data series] - 2 columns
%fs: sampling frequency in Hz
%bounds: bound in Hz for fundamental frequency range
%pSize: number of points within frequency bounds to search

%thresh: consecutive frequency values must be within thresh% of each other.
%It prohibits changes larger than thresh*mean(what_prev), where what_prev is an array of previous omega estimates. 
%Note: length of what_prev was changed from 10 to 3

%dsFactor: downsampling factor of 10 means downsample data to fs/10 Hz
%windowTime: this is the time (in seconds) duration of each window. The function does not currently allow for overlapping windows.
%R, P: harmonic and AR model orders

%A_input: transition matrix for the harmonic amplitudes (R X R). 
%If left empty, the program automatically creates an identity system state
%transition matrix (see lines 50-54 of est_harmar_emkf.m)
%sigma2w: the initial guess for s2w for the first time window only 
%omega:   the omega input allows you to manually control the fundamental frequency estimate if you already know it. 
%         If this input is zero, the program estimates omega independently at each window.
%         inputting omega as a nonzero value means that value will be interpreted as the estimate
%         Set 'omega = 0;' for automated fundamental frequency detection.

%emPeriod: this is the length of time that passes between each run of the EM s2w-update function.
%emWindow: this is the length of the window used to estimate s2w at each run of the EM algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  WINDOW PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stepTime                = windowTime;                                      %   this could be changed to allow for overlapping windows
windowSize              = fs*windowTime;                                   %   the number of samples in each fitting
emWindowSize            = fs*emWindow;                                     %   sample-length of data used in each EM estimation
stepSize                = fs*stepTime;                                     %   sample-length of steps between windows
N                       = floor(length(data)/stepSize);                    %   number of windows
emSwitch                = 1 & emPeriod;                                    %   cancels em activity if you don't want it.
Nem                     = emPeriod/windowTime;                             %   # windows between each em implementation - quick (0 if no em)
indI                    = zeros(N,1);                                      %   index where each window begins
indF                    = indI;                                            %   index where each window ends
FrqOpt                  = cell(N,1);                                       %   frequency estimates
Params                  = cell(N,1);                                       %   state parameters
Estim                   = cell(N,1);                                       %   estimated fits
ErrCov                  = cell(N,1);                                       %   error covariance outputs
NextInit                = cell(N,1);                                       %   initialization values for state estimations
EEG_Series.Raw          = data(:,2);                                       %   initialize output raw series
EEG_Series.Clean        = data(:,2);                                       %   initialize output clean EEG series
EEG_Series.BCG          = data(:,2);                                       %   initialize output BCG artifact series
data_len                = length(data);                                    %   length of data to be analyzed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FIT TO EACH WINDOW  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%wbf                     = waitbar(0,'calculating');
startI                  = 1;
endI                    = N;
for n = startI:endI
    %%   Set Indices for this Window
    disp(strcat('current window number...',num2str(n)));                   %    for diagnostics
    %pause;
    indI(n)             = 1 + (n-1)*stepSize;                              %    starting point for window n
    indF(n)             = (round(indI(n) + stepSize - 1)).*(n<endI) + (data_len).*(n==endI); %stopping point for window n
    indRange            = indI(n):1:indF(n);                               %    range of time series element to be considered in window n
    
    %%   Read in Data to Analyze in this Window
    data_fs             = [data(indRange,1) detrend(data(indRange,2))];    %   full sampling data series
    data_ds             = data_fs(1:dsFactor:end,:);                       %   downsampled series (for frequency estimation)
    
    %%   Initializations
    if n == startI                                                              
        what_prev       = [];                                              %   input parameters for estimation
        initVal         = [];                                              %   initialization parameters
        %bounds_         = bounds;                                         %   set naivee bounds
    else
        if n == startI + 1
            what_prev   = Params{n-1}.Freq*ones(3,1);                      %   use the previous frequency estimate
        else
            what_prev   = [Params{n-1}.Freq;what_prev(1:end-1)];           %   use the previous frequency estimates 
        end
        initVal         = NextInit{n-1};                                   %   initialize based on past estimates
        %bounds_        = what_prev(1)*[0.6 1.4];                          %   set bounds based on previous frequency estimates
    end
    
    %%   Run E-M Block Every emPeriod in Time
    if emSwitch && (mod(n,Nem)==startI)
        %   Update State Noise Variance s2w Estimate on block of emWindowSize seconds        
        remDataL    =  length(data) - indI(n);
        if emWindowSize <= remDataL                                        %   if there's room for a whole EM window
            nWinInEM    = floor(emWindowSize/windowSize);                  %   number of windows within the EM window
            emIndRange  = indI(n)+(0:emWindowSize-1);                      %   range of data in the window
        else                                                               
            nWinInEM    = floor(remDataL/windowSize);                      %   number of windows within the remainder of data
            emIndRange  = indI(n)+(0:nWinInEM*windowSize-1);               %   range of data in the window
            
        end
        data_fs_em      = [data(emIndRange,1) detrend(data(emIndRange,2))];%   full sampling data series
        data_ds_em      = data_fs_em(1:dsFactor:end,:);                    %   downsampled EM window
        
        %   Estimate State Noise Variance s2w
        [~,emParams,~,~,~]  = fit_harmar_emkf(data_ds_em,data_fs_em,bounds,R,P,A_input,nWinInEM,sigma2w,omega,what_prev,thresh,pSize,initVal,true, ploton);
        sigma2w         = emParams.s2w;
        s2w_iters       = emParams.s2w_iters;
        
        %   Plot convergences:
        if ploton == 1
            converge_fig = figure; set(gcf,'color','white');
        	subplot(2,1,1), plot(emParams.s2e_iters,'LineWidth',3)
            title('\^{\sigma_\epsilon^2}','FontSize',14)
            subplot(2,1,2), plot(emParams.s2w_iters,'LineWidth',3)
            title('\^{\sigma_w^2}','FontSize',14)
            pause(plot_duration);
            close(converge_fig)
        end
        clearvars emParams;
    end  
    disp(['window...' num2str(n), '...out of...',num2str(N)]);
    
    %%   Estimate Harmonic, AR and Residual Time Series Components of Data (no E-M here, so only one freq est)
    [FrqOpt{n},Params{n},Estim{n},ErrCov{n},NextInit{n}] = fit_harmar_emkf(data_ds, data_fs, bounds, R, P, A_input, 1, sigma2w, omega, what_prev,...
        thresh, pSize, initVal, false, ploton);
    EEG_Series.Raw(indRange)     = Estim{n}.Raw;
    EEG_Series.Clean(indRange)   = Estim{n}.AR;
    EEG_Series.BCG(indRange)     = Estim{n}.Harm;
    if ploton == 1
        plot_estimates(data_fs(:,1), data_fs(:,2), ...
        FrqOpt{n}, Params{n}.Freq, Params{n}.Offset, Params{n}.A, Params{n}.B, Params{n}.Alpha, Params{n}.s2e, ...
        Estim{n}.Harm, Estim{n}.ARNoise, Estim{n}.AR, Estim{n}.Noise, 4096, fs);
        pause(plot_duration); 
        close all;
    end
    
    %%   IF E-M is Running, Replace the Flat-Line with the Actual s2w Estimates
    if emSwitch && (mod(n,Nem)==startI)                    
        Params{n}.s2w_iters = s2w_iters;
    end
    
    %waitbar(n/N,wbf,['still calculating...  (' num2str(100*n/N) '%)   s2w = ' num2str(Params{n}.s2w)]);
end
%close(wbf)