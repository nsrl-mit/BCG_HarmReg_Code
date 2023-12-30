%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo for Cleaning of BCG Artifacts with DHR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add Paths
clear; clc; close all;
currDir         = [pwd, '/'];                                               %   Directory where this Script is Stored
addpath(genpath([currDir 'HarmAR']));                                       %   Code package for dynamic harmonic regression in colored noise
addpath(genpath([currDir 'Util']));                                         %   Utilies and Plotting Code Packag   

%% Set Parameters and Load Data

filename    = 'HCGSN_VEP_NoEPI.fil.bcr';                                %   Data Source to Analyze (no EPI)
data_root   = [currDir, 'data/']; % Root Folder
dataSavePath= data_root;                                                %   Path Wherein to Save Cleaned Results
load_fname  = strcat(dataSavePath, filename,'.mat');                       %   Load Data to Analyze
save_fname  = [dataSavePath filename '_clean.mat'];                     %   Filename for Results Saving
psd_fname   = [dataSavePath filename '_clean_psd.mat'];
fs          = 250;                                                      %   Sampling Frequency (Hz)

load(load_fname);                                                           %   Load Data
if ~exist('eegdata')
    eegdata     = eeg; clear eeg;                                           %   Impose Consistency of Variable Naming    
end 
if ~exist(dataSavePath,'dir')
    mkdir(dataSavePath);                                                    %   Make Directory to Save Data if it does not Exist
end                    

%% Preprocess Data
if ~exist('eegdata_filt')
H = [0 0 1 1]; f = [0 0.4/(fs/2) 0.9/(fs/2) 1]; ord = round(3*fs/0.5);      %	Filter Raw Data - Use Same Filter as fMRIB Plugin   
fwts            = firls(ord,f,H);                                           %   FIR filter weights
wo = 60/(fs/2);  Q = 35; bw = wo/Q;                                         %   Preprocessing Filter
[b,a]           = iirnotch(wo,bw);                                          %   Notch filter out power line noise  
eegdata_filt    = zeros(size(eegdata));                                     %   Initialize preprocessed EEG data stage 1(raw -> high pass filter)
eegdata_filt2   = zeros(size(eegdata));                                     %   Initialize preprocessed EEG data stage 2 (raw -> high pass filter -> notch)
for i = 1:size(eegdata,1)
    eegdata_filt(i,:) = filtfilt(fwts,1,double(eegdata(i,:)));              %   High Pass Filter EEG data  
    eegdata_filt2(i,:) = filter(b,a,eegdata_filt(i,:));                     %   Notch filter out power line noise from preprocessed EEG Data   
end
end
 
%% Clean BCG Artifacts                  
clean_chan      = [142];% , 153, 154, 163, 94, 190];                            %   List of Channels to Clean
heartbeat       = 0;                                                        %   0 (Naive HR Bounds) OR 1 (Derive HR Bounds from Events.XML)
ploton          = 0;                                                        %   Plot iterations during estimation - for diagnosis
plot_duration   = 0.5;                                                      %   Duration to keep up plots before automatically continuing (only if ploton = 1)
saveon          = 1;                                                        %   Save estimation results by channel - for comparisons across channels
[raw_eeg, clean_eeg, bcg_artifact] =  dhr_bcg_removal(eegdata_filt2', clean_chan, fs, heartbeat, ploton, dataSavePath, saveon, plot_duration);

%% Save Results
if saveon
    save(save_fname,'raw_eeg','clean_eeg','bcg_artifact','eegdata','clean_chan','fs','heartbeat','ploton','saveon');
end

%% Call Plot Diagnostics

channel_num = clean_chan(1);
logscale = 0;
plotresults = 1;

if plotresults
    plot_results_demo(save_fname)
end



