%% ADD Paths
addpath(genpath(pwd))
try
    test = itaAudio;
    clear test
catch
    disp('Indeed we need ITA-Toolbox. We shall install it, dear')
cd ..
if ispc
run([pwd '\Toolboxes\ITA-Toolbox\ita_toolbox_setup.m'])
else
    !git clone https://git.rwth-aachen.de/ita/toolbox.git ../Toolboxes/ITA-Toolbox
    run([pwd '/Toolboxes/ITA-Toolbox/ita_toolbox_setup.m'])
end

cd Calibration
end


%% Create record session
inputChannel    = 1;            %Set the channel that will receive the input 
frequencyRange  = [20 24000];   %Set the frequency range
samplingRate    = 48000;        %Set the sampling rate
microphoneCalibratorRecording = itaMSRecord('freqRange', frequencyRange,...
    'useMeasurementChain', false,'inputChannels', inputChannel,...
    'averages', 10,'fftDegree',18,'samplingRate',samplingRate);
%% Set the date
theDay = (datetime('today','format','yyyy-MM-dd'));
day = inputdlg({'day in format: dd-mm-yy:'},'Enter',[1 35],string(theDay));
%%
calibrationMicrophone = inputdlg('type yes if you want to calibrate the microphone');
%% Place the calibrator
if strcmp(calibrationMicrophone,'yes')
    disp('Place the calibrator, verify power. Use +30 dB gain. Press any key')
    pause
    % Run recording
    reference = microphoneCalibratorRecording.run;
    % iFactor is used to correct the response obtained by the microphone to Pa,
    %  making possible to know the exact amount of pressure and correcting all
    %  spectrum based on that value (Indirect calibration). That calibration is
    %  possible due to microphones' flat response
    
    % Calculate correction factor to 1 kHz
    iFactor = 1/abs(reference.freq2value(1000));
    % Verify
    signal_calibrated = (reference/2e-5)*iFactor;
    soundPressureLevel_1kHz = 20*log10(abs(signal_calibrated.freq2value(1000)));
    fprintf('Sound Pressure Level at 1000 Hz %.2f dB\n',soundPressureLevel_1kHz)
    disp('Press any key to continue')
    pause
else
    load('Current_Calibration.mat')
end

%% Measure Frequency Filters
iAverages = 3;
iFs = 48000;
iFFTDegree = 15;
iPlot = 1;
[iLoudspeakerFreqFilter,rir,rirAdjusted] = getTF(iAverages,iFs,iFFTDegree,inputChannel,iPlot);
%% Align Loudspeakers to same reproduction level
calConfig.nTolerance = 0.5;       % [dB]    Average
calConfig.nIncrement = 0.5;       % [vFS]
calConfig.nAverage = 3;           % [-]
calConfig.excitation_signal = 2;  % 1= Sweep linear || 2= pink noise || 3 = LTASS
calConfig.iChannel = inputChannel;
calConfig.nLoudspeakers = nLoudspeakers;
calConfig.level = 70;
%% playrec('reset')
[newLevelFactor,oldLevelFactor] = getLevel(iFactor,iLoudspeakerFreqFilter,calConfig);
%%
name = [day{1}];
save(['currentCalibration_' name], 'newLevelFactor','iLoudspeakerFreqFilter','iFactor','rir','rirAdjusted' )
save('current_Calibration', 'newLevelFactor','iLoudspeakerFreqFilter','iFactor','rir','rirAdjusted' )
