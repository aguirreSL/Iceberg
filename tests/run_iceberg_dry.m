% RUN_ICEBERG_DRY  Headless smoke test of iceberg_example.m
% Runs the full rendering pipeline but skips playrec playback.
% Prints a SUCCESS line + per-channel RMS so the test can be inspected
% from the MATLAB -batch output.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

try
    setToolboxes;
catch ME
    fprintf('setToolboxes warning: %s\n', ME.message);
end

calibrationPath = fullfile(projectRoot, 'src', 'calibration');
wavFilesPath    = fullfile(projectRoot, 'src', 'wavFiles');

selectedRT      = 1;
selectedAngle   = 45;
selectedLevel   = 80;
selectedSignal  = 1;   % white noise — avoids voicebox dependency for this dry run
signal = signalOptions(selectedSignal);

rtOptions = {'rt_05', 'rt_11', 'rt_00'};
selectedRT = rtOptions{selectedRT};
irFilePath = fullfile(wavFilesPath, selectedRT, 'BFormat1.wav');

[irTime, irFs] = audioread(irFilePath);
IR.time         = irTime;
IR.samplingRate = irFs;
IR.nChannels    = size(irTime, 2);
IR.nSamples     = size(irTime, 1);
IR.trackLength  = IR.nSamples / irFs;

calibrationFile = fullfile(calibrationPath, 'currentCalibration.mat');
load(calibrationFile, 'newLevelFactor', 'iFactor', 'iLoudspeakerFreqFilter');

configSetup = struct( ...
    'newLevelFactor',         newLevelFactor, ...
    'iFactor',                iFactor, ...
    'iLoudspeakerFreqFilter', iLoudspeakerFreqFilter, ...
    'wavFilesPath',           wavFilesPath);

configSetup.lsArray = [180,195,210,225,240,255,270,285,300,315,330,345, ...
                          0, 15, 30, 45, 60, 75, 90,105,120,135,150,165];
configSetup.activeLSNumbers = [1, 7, 13, 19];
iceberglAngles      = [180, 270, 0, 90];
configSetup.ls_dir  = [iceberglAngles; zeros(1,4)]';

fprintf('--- Running iceberg(...) ---\n');
iceberg_signal = iceberg(signal, IR, selectedAngle, selectedLevel, configSetup);

fprintf('\n=== SUCCESS ===\n');
fprintf('nChannels:    %d\n', iceberg_signal.nChannels);
fprintf('nSamples:     %d\n', iceberg_signal.nSamples);
fprintf('samplingRate: %d Hz\n', iceberg_signal.samplingRate);
fprintf('finite check: %d (1 = all finite)\n', all(isfinite(iceberg_signal.time(:))));
fprintf('per-channel RMS (only non-zero shown):\n');
for ch = 1:iceberg_signal.nChannels
    r = sqrt(mean(iceberg_signal.time(:, ch).^2));
    if r > 0
        fprintf('  ch %2d (angle=%3d°): RMS = %.6e\n', ch, configSetup.lsArray(ch), r);
    end
end
