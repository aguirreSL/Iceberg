% Generate MATLAB-side reference outputs for the iceberg-cpp calibration
% cross-check. For each case_X (X = 1..5):
%   - case_X/calibrate_ambisonics_ref.wav (mono output of
%     calibrate_ambisonics on input_signal.wav)
%   - case_X/calibrate_vbap_ref.wav       (4-channel output of
%     calibrate_vbap on a synthetic multichannel signal built as
%     input_signal.wav scaled by [1.0, 0.9, 0.8, 0.7])
%
% Idempotent: if the .wav refs exist and are newer than this script,
% skip the case. Editing or touching this file marks every ref stale.
%
% Run via:
%   /Applications/MATLAB_R2026a.app/bin/matlab -batch \
%     "run('Iceberg/tests/generate_calibrate_refs.m')"
%
% This script assumes iceberg-cpp-fixtures is checked out as a sibling
% of this repo:
%   git/
%     Iceberg/                <- this repo, owns every generator
%     iceberg-cpp/            <- consumes the fixtures, read-only
%     iceberg-cpp-fixtures/   <- output of this script
%
% Moved here from iceberg-cpp/tests/scripts/ on 2026-07-25: fixture
% generation is MATLAB-side work and belongs in the reference repo, so
% the fixtures directory has a single owner.

script_dir  = fileparts(mfilename('fullpath'));
projectRoot = fileparts(script_dir);
siblings    = fileparts(projectRoot);
addpath(fullfile(projectRoot, 'src', 'calibration'));
addpath(fullfile(projectRoot, 'src', 'utils'));
fixtures = fullfile(siblings, 'iceberg-cpp-fixtures');

% This script's mtime is the freshness anchor for the refs.
script_path = mfilename('fullpath');
script_info = dir([script_path '.m']);
script_mtime = script_info.datenum;

% Calibration JSON is shared across all cases.
fid = fopen([fixtures '/calibration.json']);
calib = jsondecode(fread(fid, '*char')');
fclose(fid);

ilsff = repmat(struct('freqVector', [], 'freq', []), 24, 1);
for i = 1:24
    ilsff(i).freqVector = calib.loudspeakers(i).freqVector;
    ilsff(i).freq       = calib.loudspeakers(i).freqMagnitude;
end

vbap_scales = [1.0, 0.9, 0.8, 0.7];

for ci = 1:5
    case_dir = sprintf('%s/case_%d', fixtures, ci);
    fprintf('--- case_%d ---\n', ci);

    fid = fopen([case_dir '/config.json']);
    cfg = jsondecode(fread(fid, '*char')');
    fclose(fid);

    [signal_time, fs] = audioread([case_dir '/input_signal.wav']);

    configurationSetup = struct();
    configurationSetup.iFactor                  = calib.iFactor;
    configurationSetup.newLevelFactor           = calib.newLevelFactor;
    configurationSetup.lsArray                  = cfg.lsArray;
    configurationSetup.ls_dir                   = cfg.ls_dir;
    configurationSetup.iLoudspeakerFreqFilter   = ilsff;

    % --- calibrate_ambisonics ---
    amb_ref = sprintf('%s/calibrate_ambisonics_ref.wav', case_dir);
    amb_info = dir(amb_ref);
    if ~isempty(amb_info) && amb_info.datenum > script_mtime
        fprintf('  skip ambisonics (ref is current)\n');
    else
        signal_obj = struct();
        signal_obj.time         = signal_time;
        signal_obj.samplingRate = fs;
        signal_obj.nSamples     = size(signal_time, 1);
        signal_obj.nChannels    = size(signal_time, 2);
        signal_obj.trackLength  = signal_obj.nSamples / fs;

        amb_out = calibrate_ambisonics(signal_obj, cfg.level, cfg.angle, ...
                                       configurationSetup);
        audiowrite(amb_ref, amb_out.time, fs, 'BitsPerSample', 32);
        fprintf('  wrote %s (peak=%g rms=%g)\n', amb_ref, ...
            max(abs(amb_out.time)), sqrt(mean(amb_out.time.^2)));
    end

    % --- calibrate_vbap ---
    vbap_ref = sprintf('%s/calibrate_vbap_ref.wav', case_dir);
    vbap_info = dir(vbap_ref);
    if ~isempty(vbap_info) && vbap_info.datenum > script_mtime
        fprintf('  skip vbap (ref is current)\n');
    else
        n = size(signal_time, 1);
        multich = zeros(n, 4);
        for c = 1:4
            multich(:, c) = signal_time * vbap_scales(c);
        end

        signal_obj = struct();
        signal_obj.time         = multich;
        signal_obj.samplingRate = fs;
        signal_obj.nSamples     = size(multich, 1);
        signal_obj.nChannels    = size(multich, 2);
        signal_obj.trackLength  = signal_obj.nSamples / fs;

        vbap_out = calibrate_vbap(signal_obj, cfg.level, cfg.angle, ...
                                  configurationSetup);
        audiowrite(vbap_ref, vbap_out.time, fs, 'BitsPerSample', 32);
        fprintf('  wrote %s (worst-ch peak=%g)\n', vbap_ref, ...
            max(abs(vbap_out.time(:))));
    end
end

fprintf('done.\n');
