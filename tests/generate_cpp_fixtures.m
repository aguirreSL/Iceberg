function generate_cpp_fixtures(caseNumbers)
% GENERATE_CPP_FIXTURES  Run the iceberg pipeline for selected fixture cases
% and write input/expected artifacts under:
%   <sibling of this repo>/iceberg-cpp-fixtures/case_<N>/
%
% The output directory is resolved as a sibling of this repo, matching
% generate_calibrate_refs.m and generate_iceberg_core_refs.m:
%   git/
%     Iceberg/                <- this repo, owns every generator
%     iceberg-cpp/            <- consumes the fixtures, read-only
%     iceberg-cpp-fixtures/   <- output of this script
%
% Usage:
%   generate_cpp_fixtures           % all 5 cases
%   generate_cpp_fixtures(1)        % only case 1
%   generate_cpp_fixtures([2 3 4 5])
%
% Cases:
%   1: signal=white  IR=rt_05 angle=45  level=80
%   2: signal=pink   IR=rt_05 angle=45  level=80
%   3: signal=white  IR=rt_11 angle=45  level=80
%   4: signal=white  IR=rt_05 angle= 0  level=80
%   5: signal=white  IR=rt_05 angle=90  level=80

if nargin < 1 || isempty(caseNumbers)
    caseNumbers = 1:5;
end

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
try
    setToolboxes;
catch ME
    fprintf('setToolboxes warning: %s\n', ME.message);
end

calibrationPath = fullfile(projectRoot, 'src', 'calibration');
wavFilesPath    = fullfile(projectRoot, 'src', 'wavFiles');
outRoot         = fullfile(fileparts(projectRoot), 'iceberg-cpp-fixtures');
if ~exist(outRoot, 'dir'); mkdir(outRoot); end

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
icebergAngles      = [180, 270, 0, 90];
configSetup.ls_dir = [icebergAngles; zeros(1,4)]';

cases = { ...
    struct('id',1,'signal_id',1,'ir_name','rt_05','angle',45,'level',80), ...
    struct('id',2,'signal_id',2,'ir_name','rt_05','angle',45,'level',80), ...
    struct('id',3,'signal_id',1,'ir_name','rt_11','angle',45,'level',80), ...
    struct('id',4,'signal_id',1,'ir_name','rt_05','angle', 0,'level',80), ...
    struct('id',5,'signal_id',1,'ir_name','rt_05','angle',90,'level',80) };

for k = caseNumbers(:)'
    c = cases{k};
    fprintf('\n=== CASE %d: signal=%d ir=%s angle=%d level=%d ===\n', ...
        c.id, c.signal_id, c.ir_name, c.angle, c.level);

    rng(42);
    signal = signalOptions(c.signal_id);

    irFilePath = fullfile(wavFilesPath, c.ir_name, 'BFormat1.Wav');
    [irTime, irFs] = audioread(irFilePath);
    if irFs ~= signal.samplingRate
        error('Sample rate mismatch: signal=%d Hz, IR(%s)=%d Hz. Resample required.', ...
            signal.samplingRate, c.ir_name, irFs);
    end
    IR = struct('time', irTime, 'samplingRate', irFs, ...
                'nChannels', size(irTime,2), 'nSamples', size(irTime,1));
    IR.trackLength = IR.nSamples / IR.samplingRate;

    iceberg_signal = iceberg(signal, IR, c.angle, c.level, configSetup);

    outDir = fullfile(outRoot, sprintf('case_%d', c.id));
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    % 32-bit float WAV — output magnitudes are voltage-scaled, can exceed [-1,1]
    audiowrite(fullfile(outDir,'input_signal.wav'),    signal.time,         signal.samplingRate,         'BitsPerSample', 32);
    audiowrite(fullfile(outDir,'input_ir.wav'),        IR.time,             IR.samplingRate,             'BitsPerSample', 32);
    audiowrite(fullfile(outDir,'expected_output.wav'), iceberg_signal.time, iceberg_signal.samplingRate, 'BitsPerSample', 32);

    nCh = iceberg_signal.nChannels;
    rms_per_ch  = zeros(1, nCh);
    peak_per_ch = zeros(1, nCh);
    for ch = 1:nCh
        x = iceberg_signal.time(:, ch);
        rms_per_ch(ch)  = sqrt(mean(x.^2));
        peak_per_ch(ch) = max(abs(x));
    end

    cfg = struct();
    cfg.case_id          = c.id;
    cfg.signal_id        = c.signal_id;
    cfg.ir_name          = c.ir_name;
    cfg.angle            = c.angle;
    cfg.level            = c.level;
    cfg.lsArray          = configSetup.lsArray;
    cfg.activeLSNumbers  = configSetup.activeLSNumbers;
    cfg.ls_dir           = configSetup.ls_dir;
    cfg.samplingRate     = iceberg_signal.samplingRate;
    cfg.signal_nSamples  = signal.nSamples;
    cfg.ir_nSamples      = IR.nSamples;
    cfg.ir_nChannels     = IR.nChannels;
    cfg.rng_seed         = 42;
    fid = fopen(fullfile(outDir,'config.json'), 'w');
    fwrite(fid, jsonencode(cfg));
    fclose(fid);

    metrics = struct();
    metrics.case_id          = c.id;
    metrics.nSamples         = iceberg_signal.nSamples;
    metrics.nChannels        = iceberg_signal.nChannels;
    metrics.samplingRate     = iceberg_signal.samplingRate;
    metrics.rms_per_channel  = rms_per_ch;
    metrics.peak_per_channel = peak_per_ch;
    fid = fopen(fullfile(outDir,'expected_metrics.json'), 'w');
    fwrite(fid, jsonencode(metrics));
    fclose(fid);

    fprintf('Wrote %s/{input_signal.wav,input_ir.wav,expected_output.wav,config.json,expected_metrics.json}\n', outDir);
    fprintf('Active-channel metrics (ch idx | angle | RMS | peak):\n');
    for ch = configSetup.activeLSNumbers
        fprintf('  ch %2d (angle=%3d°): RMS=%.6e  peak=%.6e\n', ...
            ch, configSetup.lsArray(ch), rms_per_ch(ch), peak_per_ch(ch));
    end
end
end
