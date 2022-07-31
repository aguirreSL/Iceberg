%% Add Paths/ITA

addpath(genpath(pwd))
setToolboxes
try
    test = itaAudio;
    clear test
catch
    disp('Indeed we need ITA-Toolbox. We shall install it, dear')
    if ispc
    run([pwd '\Toolboxes\ITA-Toolbox\ita_toolbox_setup.m'])
    elseif ismac
    run([pwd '\Toolboxes\ITA-Toolbox\ita_toolbox_setup.m'])
    end
end

ccx % Clear all, clean all, close all :D
calibrationPath = pwd;
if ispc
    load([calibrationPath '\Current_Calibration.mat'], 'new_Level_Factor','iFactor','iLoudspeakerFreqFilter');
elseif ismac
    load([calibrationPath '/Current_Calibration.mat'], 'new_Level_Factor','iFactor','iLoudspeakerFreqFilter');
end
    
%% Calibration example file
% % load('playpen_calibration_11-Sep-2018.mat', 'full_calibration')

%These are values to make it work. You should run a calibration session to
%your specific array
configurationSetup.new_Level_Factor = new_Level_Factor;
configurationSetup.iFactor = iFactor;
configurationSetup.iLoudspeakerFreqFilter = iLoudspeakerFreqFilter;

%% Array Specific Settings. 
% Set here the channels relative to the angles that follow the input channel number of the sound card
% This script assumes the Glasgow (River Clyde room) with 24 LS and clockwise configuration:
% (LS 1 is the backward one (180 degrees) / LS 13 is the frontal (0 degrees))
% (I'm assuming the same configuration to UoN)

Channel_1 =  180; %The back loudspeaker is connected to output number 1 of the sound card. 
Channel_7 =  270; %The left loudspeaker is connected to output number 7 of the sound card. 
Channel_13 = 000; %The front loudspeaker is connected to output number 13 of the sound card. 
Channel_19 = 090; %The right loudspeaker is connected to output number 19 of the sound card. 

TotalNumberOfLoudspeakers = 24;

spaceBetweenLS = 360/TotalNumberOfLoudspeakers;
% To be used in 4 LS Hybrid mode we need to specify the LS index numbers
configurationSetup.ls_dir = [Channel_1 Channel_7 Channel_13 Channel_19; 0 0 0 0]'; 
% Clockwise configuration (I'm assuming the same configuration to UoN)
% Glasgow_array = [180:15:345 0:15:165]'; 
%Right LS is 90 Degrees / Left is 270
configurationSetup.LSArray = [mod(Channel_1:spaceBetweenLS:360,360) spaceBetweenLS:spaceBetweenLS:(Channel_1-spaceBetweenLS)]'; %Clockwise configuration
%Right LS is 270 Degrees / Left is 90 counter = clockwise configuration
% configurationSetup.LSArray = [mod(Channel_1:-spaceBetweenLS:0,360) (Channel_1-spaceBetweenLS):-spaceBetweenLS:spaceBetweenLS]'; 
sampleFrequency = 44100; %Upt o 96k with Motu/ Maybe ferrofish as well
lengthSignalSeconds = 3; %Set the desired Lenght here

fftDegree = log2(lengthSignalSeconds*44100); %From seconds to fft degree

frequencyLimits = [20 20000]; % probably the tanoy speakers are fine from 50/60 Hz with the inverted filter. I've added to my todo list []check FRF from Glasgow LS

%% Signals (I've created some options, you can define/import any signal) %%% 
% signal_1 White Noise
signal_1 = ita_generate('whitenoise',1,sampleFrequency,fftDegree);
% signal_2 Pink Noise
signal_2 = ita_generate('pinknoise',1,sampleFrequency,fftDegree);
% signal_3 LTASS - long-term average speech spectrum
if ismac
    file = ita_read(strrep('.\wavFiles\HINT_Danish_Noise_Female1_Inf_new.wav','\','/'));             %importing the FEMALE speech signal
else
    file = ita_read('.\wavFiles\HINT_Danish_Noise_Female1_Inf_new.wav');             %importing the FEMALE speech signal
end
    [vs,~]=vadsohn(file.time,file.samplingRate,'p');                                %Voice Active Detector
    zeroing = [(vs'.*file.time(1:length(vs))'),0];                                  %Mute silent parts  
    file = itaAudio([nonzeros(zeroing);0],file.samplingRate,'time');                %ITA Audio
    smoothSignal        = ita_smooth_frequency(file);                             %Smoth signal (Frequency shape)   
    whiteNoiseSignal    = ita_generate('noise',1,file.samplingRate,file.fftDegree); %Generate white-noise signal/ (is pink more appropriated?)
    LTASS            = ita_multiply_spk(whiteNoiseSignal,smoothSignal);             %Shaping random signal to speech
    LTASS            = ita_time_crop(LTASS,[0 lengthSignalSeconds],'time');         %Defined time (it should be smaller than the imported speech audio)
    LTASS            = ita_normalize_dat(LTASS);                                    %Normalize time domain
    signal_3 = ita_resample(LTASS,sampleFrequency);                                 %Resampling if needed
% signal_4 Log sweep (if more energy in lower frequency is required)
signal_4 = ita_generate('ccxsweep',frequencyLimits,sampleFrequency,fftDegree);
% signal_5 Linear sweep
signal_5 = ita_generate('swenlinsweep',frequencyLimits,0.0,sampleFrequency,fftDegree);
% signal_6 ISTS International Speech Test Signal (repeated if statement to improve readability)
if ismac
    signal_6 = ita_read([pwd strrep('\wavFiles\ISTS-V1.0_60s_16bit.wav','\','/')]);
else
    signal_6 = ita_read([pwd '\wavFiles\ISTS-V1.0_60s_16bit.wav']);
end
%% Select the reverberation time, presentation angle, presentation level, and the signal. 
Play_At_This_RT = 3; %1 = 0.5 s 
                     %2 = 1.1 s
                     %3 = 0.0 s
Play_From_This_Angle    = 45; %0:5:355
Play_At_This_Level      = 85; %dB SPL
Play_This_Signal        = signal_3; 

%% Creating the signal

[Signal_To_Run] = setAuralizationHybrid_4LS_2021(Play_This_Signal,Play_At_This_Level,Play_From_This_Angle,Play_At_This_RT,configurationSetup)
  






