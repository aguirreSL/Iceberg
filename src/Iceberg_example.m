%% Welcome to the Iceberg Auralization Method
%This is an example on how to create an auralized file using a combination
%of VBAP and Ambisonics. The IR is divided using Room Acoustics Center Time.
%Three examples of RT were made available (0, 0.5, and 1.1 s)
addpath(genpath(pwd))
setToolboxes
ccx % Clear all, clean all, close all :D
%% Calibration example file
calibrationPath = fullfile(pwd, 'Calibration');
toolboxesPath = fullfile(pwd, 'Toolboxes');
wavFilesPath = fullfile(pwd, 'wavFiles');
currentDir = pwd;
sep = '/';
if ispc
    sep = '\';  
end
load([currentDir sep 'Calibration' sep 'Current_Calibration.mat'], 'new_Level_Factor','iFactor','iLoudspeakerFreqFilter');
%These are placeholder values. 
%You should run a calibration session to your specific setup

configurationSetup.new_Level_Factor = new_Level_Factor;
configurationSetup.iFactor = iFactor;
configurationSetup.iLoudspeakerFreqFilter = iLoudspeakerFreqFilter;
configurationSetup.wavFilesPath = wavFilesPath;

%% Array Specific Settings. 
% Set here the channels relative to the angles that follow the input channel number of the sound card
% This script assumes the Glasgow (River Clyde room) with 24 LS and clockwise configuration:
% (LS 1 is the backward one (180 degrees) / LS 13 is the frontal (0 degrees))

TotalNumberOfLoudspeakers = 24;
spaceBetweenLS = 360/TotalNumberOfLoudspeakers;

Channel_1 =  180; %The back loudspeaker is connected to output number 1 of the sound card. 
Channel_7 =  270; %The left loudspeaker is connected to output number 7 of the sound card. 
Channel_13 = 000; %The front loudspeaker is connected to output number 13 of the sound card. 
Channel_19 = 090; %The right loudspeaker is connected to output number 19 of the sound card. 

% To be used in 4 LS Hybrid mode we need to specify the LS index numbers
configurationSetup.ls_dir = [Channel_1 Channel_7 Channel_13 Channel_19; 0 0 0 0]'; 
%Right LS is 90 Degrees / Left is 270 %Clockwise configuration
configurationSetup.LSArray = [mod(Channel_1:spaceBetweenLS:360,360) spaceBetweenLS:spaceBetweenLS:(Channel_1-spaceBetweenLS)]'; 
%%[option] Right LS is 270 Degrees / Left is 90 counter = clockwise configuration
% configurationSetup.LSArray = [mod(Channel_1:-spaceBetweenLS:0,360) (Channel_1-spaceBetweenLS):-spaceBetweenLS:spaceBetweenLS]';

%% Select the reverberation time, Presentation angle, Presentation Level, and the Signal. 
Selected_RT      = 1;  % 1 = 0.5 s | %2 = 1.1 s | %3 = 0.0 s
Selected_Angle   = 45; % 0:5:355
Selected_Level   = 80; % dB SPL
Selected_Signal  = 7;  % options 1 to 7 described in the 'signalOptions' function
selectedSignal = signalOptions(Selected_Signal);

%% Load the RIR (Odeon Simulations, can be any 1st order ambisonics)
if Selected_RT == 1
    selectRT = 'rt_05';
elseif Selected_RT == 2
    selectRT = 'rt_11';
elseif Selected_RT == 3
    selectRT = 'rt_00';
end

IR = ita_read([wavFilesPath sep selectRT sep 'BFormat1.Wav']);   %Load Ambisonics IR

[DSER,LR] = iceberg_core(IR);

[iceberg_signal, VBAP_Part, Amb_Part] = iceberg_merge(selectedSignal, DSER, LR, ...
           Selected_Level, Selected_Angle, configurationSetup);

%% Play the resulting file
%Set you soundcard in ita_preferences first
playPlain(iceberg_signal)

%% Play function
function playPlain(signal)
[~, playDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('playDeviceID'));
[~, recDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('recDeviceID'));
playDeviceID = playDeviceInfo.deviceID;
recDeviceID = recDeviceInfo.deviceID;
iFs = signal.samplingRate;                %Sample Frequency
signal_run = signal;
nCh = 24;
if playrec('isInitialised')==1
    playrec('reset')
    playrec('init', iFs, playDeviceID, recDeviceID)
    playrec('play',signal_run.time,1:nCh);
else
    playrec('init', iFs, playDeviceID, recDeviceID)
    playrec('play',signal_run.time,1:nCh);
end
end