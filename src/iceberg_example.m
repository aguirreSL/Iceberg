%% ICEBERG AURALIZATION METHOD: VBAP-AMBISONICS HYBRID APPROACH
% This script implements a hybrid spatial audio rendering technique combining 
% Vector Base Amplitude Panning (VBAP) for direct sound and early reflections
% with 1st Order Ambisonics for late reverberation.
%
% Key Features:
%   - Uses Room Acoustics Center Time to divide impulse responses
%   - Supports three reverberation times: 0.0s, 0.5s, and 1.1s
%   - Configurable presentation angle, level, and source signal
%
% Dependencies:
%   - ITA-Toolbox
%   - polarch HOA, VBAP and SAP-Voicebox toolboxes
%   - Supporting custom functions
%
% System Requirements:
%   - MATLAB R2018b or later
%   - ASIO-compatible sound card for multichannel playback
%
% Updated: 09/06/2025

%% Initialization
addpath(genpath(pwd));      % Add all subdirectories to path
setToolboxes;               % Initialize required toolboxes
ccx;                        % Clear workspace, command window, and close figures
calibrationPath = fullfile(pwd, 'Calibration');
wavFilesPath    = fullfile(pwd, 'wavFiles');

%% AURALIZATION PARAMETERS
selectedRT      = 1;     % 1:0.5s, 2:1.1s, 3:0.0s
selectedAngle   = 45;    % Presentation angle [0-355°]
selectedLevel   = 80;    % Output level [dB SPL]
selectedSignal  = 3;     % Source signal option (see signalOptions.m)
signal = signalOptions(selectedSignal);

%% LOAD IMPULSE RESPONSE
rtOptions = {'rt_05', 'rt_11', 'rt_00'};
selectedRT = rtOptions{selectedRT};
irFilePath = fullfile(wavFilesPath, selectedRT, 'BFormat1.wav');
IR = ita_read(irFilePath);
% Load calibration data (placeholder values)
%You should run a calibration session to your specific setup
calibrationFile = fullfile(calibrationPath, 'Current_Calibration.mat');
load(calibrationFile, 'new_Level_Factor', 'iFactor', 'iLoudspeakerFreqFilter');

configSetup = struct(...
    'new_Level_Factor', new_Level_Factor,...
    'iFactor', iFactor,...
    'iLoudspeakerFreqFilter', iLoudspeakerFreqFilter,...
    'wavFilesPath', wavFilesPath);

clear new_Level_Factor iFactor iLoudspeakerFreqFilter

%% ARRAY CONFIGURATION SETTINGS
% Defines the channel-to-angle mapping for the loudspeaker array.
% Current setup: Glasgow (River Clyde room) 24-loudspeaker circular array
% - Array is arranged in clockwise configuration
% - Channel 1:  180° (back)  |  - Channel  7:  270° (left)
% - Channel 13:   0° (front) |  - Channel 19:   90°  (right)
numberOfSpeakers = 24;
configSetup.lsArray = [180,195,210,225,240,255,270,285,300,315,330,345,...
                            0,15,30,45,60,75,90,105,120,135,150,165];
% Define 4 LS to be used by Iceberg
configSetup.activeLSNumbers = [1, 7, 13, 19];
iceberglAngles           = [180, 270, 0, 90]; % [Back, Left, Front, Right]
configSetup.ls_dir       = [iceberglAngles; zeros(1,4)]';

%% PROCESS SIGNAL
iceberg_signal = iceberg(signal, IR, selectedAngle, configSetup);

%% PLAYBACK
playPlain(iceberg_signal,numberOfSpeakers)% Output to configured sound card

%% PLAYBACK FUNCTION
function playPlain(signal_run,nCh)
[~, playDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('playDeviceID'));
[~, recDeviceInfo]  = ita_portaudio_deviceID2string(ita_preferences('recDeviceID'));
if playrec('isInitialised')==1
    playrec('reset')
end
    playrec('init', signal_run.samplingRate, playDeviceInfo.deviceID, recDeviceInfo.deviceID)
    playrec('play', signal_run.time,1:nCh);
end