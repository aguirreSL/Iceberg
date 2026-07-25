function export_calibration()
% EXPORT_CALIBRATION  Dump src/calibration/currentCalibration.mat to JSON for
% the C++ port. Writes to:
%   <sibling of this repo>/iceberg-cpp-fixtures/calibration.json
%
% The output directory is resolved as a sibling of this repo, matching the
% three generators under tests/:
%   git/
%     Iceberg/                <- this repo, owns every generator
%     iceberg-cpp/            <- consumes the fixtures, read-only
%     iceberg-cpp-fixtures/   <- output of this script
%
% Run from the repo root:
%   matlab -batch "addpath(genpath(pwd)); export_calibration"

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
try
    setToolboxes;
catch ME
    fprintf('setToolboxes warning: %s\n', ME.message);
end

calibrationFile = fullfile(projectRoot, 'src', 'calibration', 'currentCalibration.mat');
S = load(calibrationFile, 'newLevelFactor', 'iFactor', 'iLoudspeakerFreqFilter');

iLoudspeakerFreqFilter = S.iLoudspeakerFreqFilter;
newLevelFactor         = S.newLevelFactor;
iFactor                = S.iFactor;

nLS = length(iLoudspeakerFreqFilter);
if nLS ~= 24
    warning('Expected 24 loudspeakers, got %d', nLS);
end

loudspeakers = cell(1, nLS);
for i = 1:nLS
    filt = iLoudspeakerFreqFilter(i);
    fv   = double(filt.freqVector);
    fm   = double(abs(filt.freq));
    loudspeakers{i} = struct( ...
        'index',         i, ...
        'freqVector',    fv(:).', ...
        'freqMagnitude', fm(:).');
end

out = struct();
out.iFactor        = double(iFactor);
out.newLevelFactor = double(newLevelFactor(:).');
out.loudspeakers   = loudspeakers;

outDir  = fullfile(fileparts(projectRoot), 'iceberg-cpp-fixtures');
if ~exist(outDir, 'dir'); mkdir(outDir); end
outPath = fullfile(outDir, 'calibration.json');

fid = fopen(outPath, 'w');
if fid == -1
    error('Cannot open output for writing: %s', outPath);
end
fwrite(fid, jsonencode(out));
fclose(fid);

fprintf('Wrote %s\n', outPath);
fprintf('  iFactor                 = %.6g\n', out.iFactor);
fprintf('  newLevelFactor length   = %d\n', numel(out.newLevelFactor));
fprintf('  loudspeakers count      = %d\n', nLS);
fprintf('  ls(1).freqVector length = %d\n', numel(loudspeakers{1}.freqVector));
fprintf('  ls(1).freqMagnitude len = %d\n', numel(loudspeakers{1}.freqMagnitude));
end
