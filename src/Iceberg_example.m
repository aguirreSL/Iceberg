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
if ispc
    load([currentDir '\Calibration\Current_Calibration.mat'], 'new_Level_Factor','iFactor','iLoudspeakerFreqFilter');
elseif ismac
    load([currentDir '/Calibration/Current_Calibration.mat'], 'new_Level_Factor','iFactor','iLoudspeakerFreqFilter');
end

%These are values to make it work. You should run a calibration session to
%your specific array
configurationSetup.new_Level_Factor = new_Level_Factor;
configurationSetup.iFactor = iFactor;
configurationSetup.iLoudspeakerFreqFilter = iLoudspeakerFreqFilter;
configurationSetup.wavFilesPath = wavFilesPath;

%% Array Specific Settings. 
% Set here the channels relative to the angles that follow the input channel number of the sound card
% This script assumes the Glasgow (River Clyde room) with 24 LS and clockwise configuration:
% (LS 1 is the backward one (180 degrees) / LS 13 is the frontal (0 degrees))

Channel_1 =  180; %The back loudspeaker is connected to output number 1 of the sound card. 
Channel_7 =  270; %The left loudspeaker is connected to output number 7 of the sound card. 
Channel_13 = 000; %The front loudspeaker is connected to output number 13 of the sound card. 
Channel_19 = 090; %The right loudspeaker is connected to output number 19 of the sound card. 

TotalNumberOfLoudspeakers = 24;

spaceBetweenLS = 360/TotalNumberOfLoudspeakers;
% To be used in 4 LS Hybrid mode we need to specify the LS index numbers
configurationSetup.ls_dir = [Channel_1 Channel_7 Channel_13 Channel_19; 0 0 0 0]'; 
% Glasgow_array = [180:15:345 0:15:165]'; 
%Right LS is 90 Degrees / Left is 270
configurationSetup.LSArray = [mod(Channel_1:spaceBetweenLS:360,360) spaceBetweenLS:spaceBetweenLS:(Channel_1-spaceBetweenLS)]'; %Clockwise configuration
%Right LS is 270 Degrees / Left is 90 counter = clockwise configuration
% configurationSetup.LSArray = [mod(Channel_1:-spaceBetweenLS:0,360) (Channel_1-spaceBetweenLS):-spaceBetweenLS:spaceBetweenLS]'; 


%% Select the reverberation time, Presentation angle, Presentation Level, and the Signal. 
Play_At_This_RT = 3; %1 = 0.5 s 
                     %2 = 1.1 s
                     %3 = 0.0 s
Play_From_This_Angle    = 45; %0:5:355
Play_At_This_Level      = 85; %dB SPL
Play_This_Signal        = 7; %options 1 to 7 described in the next section

%% Signals (I've created some options, you can define/import any signal) %%% 
sampleFrequency = 44100; %Up to 96k with Motu/ Maybe ferrofish as well
lengthSignalSeconds = 3; %Set the desired Lenght here
fftDegree = log2(lengthSignalSeconds*sampleFrequency); %From seconds to fft degree
frequencyLimits = [20 20000]; % probably the tanoy speakers are fine from 50/60 Hz with the inverted filter. I've added to my todo list []check FRF from Glasgow LS

switch Play_This_Signal
    case 1
        % signal_1 White Noise
        selectedSignal = ita_generate('whitenoise',1,sampleFrequency,fftDegree);        
    case 2
        % signal_2 Pink Noise
        selectedSignal = ita_generate('pinknoise',1,sampleFrequency,fftDegree);
    case 3
        %signal_3 LTASS - long-term average speech spectrum
        %     file = ita_read('HINT_Danish_Noise_Male1_Inf_new.wav');                       %importing the MALE speech signal 
            file = ita_read('.\wavFiles\HINT_Danish_Noise_Female1_Inf_new.wav');             %importing the FEMALE speech signal
            [vs,~]=vadsohn(file.time,file.samplingRate,'p');                                %Voice Active Detector
            zeroing = [(vs'.*file.time(1:length(vs))'),0];                                  %Mute silent parts  
            file = itaAudio([nonzeros(zeroing);0],file.samplingRate,'time');                %ITA Audio
            smoothSignal        = ita_smooth_frequency(file);                             %Smoth signal (Frequency shape)   
            whiteNoiseSignal    = ita_generate('noise',1,file.samplingRate,file.fftDegree); %Generate white-noise signal/ (is pink more appropriated?)
            LTASS            = ita_multiply_spk(whiteNoiseSignal,smoothSignal);             %Shaping random signal to speech
            LTASS            = ita_time_crop(LTASS,[0 lengthSignalSeconds],'time');         %Defined time (it should be smaller than the imported speech audio)
            LTASS            = ita_normalize_dat(LTASS);                                    %Normalize time domain
            selectedSignal = ita_resample(LTASS,sampleFrequency);                                 %Resampling if needed
    case 4
        % signal_4 Log sweep
        selectedSignal = ita_generate('ccxsweep',frequencyLimits,sampleFrequency,fftDegree);
    case 5
        % signal_5 Linear sweep
        selectedSignal = ita_generate('swenlinsweep',frequencyLimits,0.0,sampleFrequency,fftDegree);
    case 6
        % signal_6 ISTS International Speech Test Signal
        selectedSignal = ita_read([pwd '\wavFiles\ISTS-V1.0_60s_16bit.wav']);
    case 7
        %Pure tone 400 Hz 
        selectedSignal = ita_generate('sine',1,400,sampleFrequency,fftDegree);
end

%% Creating the signal

[Signal_To_Run] = Iceberg(selectedSignal,Play_At_This_Level,Play_From_This_Angle,Play_At_This_RT,configurationSetup)
  
%% Play the resulting file

playPlain(Signal_To_Run)

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





