function iceberg_signal = iceberg(signal, IR, Selected_Angle, level, configSetup)
% ICEBERG Auralization Method (VBAP + Ambisonics)
% Native Implementation
% Combines VBAP for direct sound and early reflections with Ambisonics for
% late reverberation in spatial audio rendering.
%
% Inputs:
%   signal          - Input audio Struct (with .time, .samplingRate)
%   IR              - Impulse responses Struct
%   Selected_Angle  - Source presentation angle (degrees)
%   level           - Target SPL in dB (or 'n' to bypass level scaling)
%   configSetup     - Configuration setup structure
%
% Output:
%   iceberg_signal  - Processed spatial audio Struct

[dser, lr] = iceberg_core(IR);

vbapDser = iceberg_set_vbap(signal, dser, Selected_Angle, level, configSetup);

ambisonicsLr = iceberg_set_amb(signal, lr, Selected_Angle, level, configSetup);

% Pre-populate empty audio structs for the active channels
vbapAudio.time = zeros(size(vbapDser.time, 1), max(configSetup.activeLSNumbers));
vbapAudio.samplingRate = signal.samplingRate;

ambAudio.time = zeros(size(ambisonicsLr.time, 1), max(configSetup.activeLSNumbers));
ambAudio.samplingRate = signal.samplingRate;

activeLsNumbers = configSetup.activeLSNumbers;

% Distribute signals to appropriate loudspeaker channels
for i = 1:length(activeLsNumbers)
    channel = activeLsNumbers(i);
    vbapAudio.time(:, channel) = vbapDser.time(:, i);
    
    % Only write to Ambisonics if it output enough channels
    if i <= size(ambisonicsLr.time, 2)
        ambAudio.time(:, channel)  = ambisonicsLr.time(:, i);
    end
end

% update native structures dimensions
vbapAudio.nChannels = size(vbapAudio.time, 2);
vbapAudio.nSamples = size(vbapAudio.time, 1);
ambAudio.nChannels = size(ambAudio.time, 2);
ambAudio.nSamples = size(ambAudio.time, 1);

% Merge together
iceberg_signal = native_add(vbapAudio, ambAudio);

% Pad array with zeros up to the full array layout
targetChannels = length(configSetup.lsArray);
currentChannels = iceberg_signal.nChannels;

if targetChannels > currentChannels
    iceberg_signal.time(:, currentChannels+1:targetChannels) = ...
        zeros(iceberg_signal.nSamples, targetChannels - currentChannels);
    iceberg_signal.nChannels = size(iceberg_signal.time, 2);
end

end