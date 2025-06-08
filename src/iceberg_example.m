%% ICBERG AURALIZATION METHOD: VBAP-AMBISONICS HYBRID APPROACH
% This script demonstrates how to create an auralized audio file using a hybrid 
% approach combining Vector Base Amplitude Panning (VBAP) and Ambisonics.
%
% Key Features:
%   - Uses Room Acoustics Center Time to divide impulse responses
%   - Supports three reverberation times: 0.0s, 0.5s, and 1.1s
%   - Configurable presentation angle, level, and source signal
%
% Required Toolboxes: ITA-Toolbox, polarch hoa, polarch vbap, and supporting functions
% Updated: 2025-06-08
%% Initialization
addpath(genpath(pwd));      % Add all subdirectories to path
setToolboxes;               % Initialize required toolboxes
ccx;                        % Clear workspace, command window, and close figures

%% Path Configuration
calibrationPath = fullfile(pwd, 'Calibration');
toolboxesPath = fullfile(pwd, 'Toolboxes');
wavFilesPath = fullfile(pwd, 'wavFiles');

% Load calibration data (placeholder values)
%You should run a calibration session to your specific setup
calibrationFile = fullfile(calibrationPath, 'Current_Calibration.mat');
load(calibrationFile, 'new_Level_Factor', 'iFactor', 'iLoudspeakerFreqFilter');

configurationSetup = struct(...
    'new_Level_Factor', new_Level_Factor,...
    'iFactor', iFactor,...
    'iLoudspeakerFreqFilter', iLoudspeakerFreqFilter,...
    'wavFilesPath', wavFilesPath);

%% ARRAY CONFIGURATION SETTINGS
% Defines the channel-to-angle mapping for the loudspeaker array.
% Current setup: Glasgow (River Clyde room) 24-loudspeaker circular array
% - Array is arranged in clockwise configuration
% - Channel 1:  180° (back)
% - Channel 7:  270° (left)
% - Channel 13: 0°   (front)
% - Channel 19: 90°  (right)
% Note: Channel numbers correspond to sound card input channels

TotalNumberOfLoudspeakers = 24;
angularSpacing = 360/TotalNumberOfLoudspeakers;

% Define cardinal loudspeaker positions
cardinalChannels = [1, 7, 13, 19];
cardinalAngles = [180, 270, 0, 90];  % [Back, Left, Front, Right]
configurationSetup.ls_dir = [cardinalAngles; zeros(1,4)]';

% Generate loudspeaker positions starting from Channel 1 (180°)
% Clockwise arrangement: 180° ? 195° ? ... ? 0° ? 15° ? ... ? 165°
lsAnglesPart1 = mod(cardinalAngles(1):angularSpacing:360, 360);  % 180° to 360° (0°)
lsAnglesPart2 = angularSpacing:angularSpacing:165;  % 15° to 165°
configurationSetup.LSArray = [lsAnglesPart1, lsAnglesPart2]';
%% AURALIZATION PARAMETERS
% Select processing options:
Selected_RT      = 1;     % 1:0.5s, 2:1.1s, 3:0.0s
Selected_Angle   = 45;    % Presentation angle [0-355°]
Selected_Level   = 80;    % Output level [dB SPL]
Selected_Signal  = 7;     % Source signal option (see signalOptions.m)

%% LOAD IMPULSE RESPONSE
% Select corresponding reverberation time dataset
rtOptions = {'rt_05', 'rt_11', 'rt_00'};
selectedRT = rtOptions{Selected_RT};

% Load B-Format Ambisonics impulse response
irFilePath = fullfile(wavFilesPath, selectedRT, 'BFormat1.wav');
IR = ita_read(irFilePath);
       
%% PROCESS SIGNAL
% 1. Perform core Iceberg processing
[DSER, LR] = iceberg_core(IR);

% 2. Get selected audio signal
selectedSignal = signalOptions(Selected_Signal);

% 3. Generate hybrid VBAP/Ambisonics output
[iceberg_signal, VBAP_Part, Amb_Part] = iceberg_merge(...
    selectedSignal, DSER, LR, Selected_Level, Selected_Angle, configurationSetup);       

%% PLAYBACK
playPlain(iceberg_signal,TotalNumberOfLoudspeakers)% Output to configured sound card

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