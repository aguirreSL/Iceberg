% Generate MATLAB-side references for iceberg_core on all 5 iceberg-cpp
% fixtures. Outputs:
%   case_X/iceberg_core_dser_ref.wav    (mono DSER)
%   case_X/iceberg_core_ir_late_ref.wav (4-ch IR_Late)
%
% Idempotent: skips when refs are newer than this script. Editing or
% touching this file marks every ref stale.
%
% Run via:
%   matlab -batch "run('Iceberg/tests/generate_iceberg_core_refs.m')"
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
addpath(fullfile(projectRoot, 'src', 'core'));
addpath(fullfile(projectRoot, 'src', 'utils'));
fixtures = fullfile(siblings, 'iceberg-cpp-fixtures');

script_path = mfilename('fullpath');
script_info = dir([script_path '.m']);
script_mtime = script_info.datenum;

for ci = 1:5
    case_dir = sprintf('%s/case_%d', fixtures, ci);
    fprintf('--- case_%d ---\n', ci);

    [ir_time, fs] = audioread([case_dir '/input_ir.wav']);

    ir_struct = struct();
    ir_struct.time         = ir_time;
    ir_struct.samplingRate = fs;
    ir_struct.nSamples     = size(ir_time, 1);
    ir_struct.nChannels    = size(ir_time, 2);
    ir_struct.trackLength  = ir_struct.nSamples / fs;

    dser_ref    = sprintf('%s/iceberg_core_dser_ref.wav', case_dir);
    ir_late_ref = sprintf('%s/iceberg_core_ir_late_ref.wav', case_dir);

    dser_info    = dir(dser_ref);
    ir_late_info = dir(ir_late_ref);
    if ~isempty(dser_info) && dser_info.datenum > script_mtime && ...
       ~isempty(ir_late_info) && ir_late_info.datenum > script_mtime
        fprintf('  skip (refs current)\n');
        continue;
    end

    [DSER, IR_Late] = iceberg_core(ir_struct);

    audiowrite(dser_ref,    DSER.time,    fs, 'BitsPerSample', 32);
    audiowrite(ir_late_ref, IR_Late.time, fs, 'BitsPerSample', 32);

    fprintf('  DSER   : %d x %d  peak=%.6g rms=%.6g\n', ...
        size(DSER.time, 1), size(DSER.time, 2), ...
        max(abs(DSER.time(:))), sqrt(mean(DSER.time(:).^2)));
    fprintf('  IR_Late: %d x %d  peak=%.6g rms=%.6g\n', ...
        size(IR_Late.time, 1), size(IR_Late.time, 2), ...
        max(abs(IR_Late.time(:))), sqrt(mean(IR_Late.time(:).^2)));
end

fprintf('done.\n');
