function signal_ambisonics = iceberg_set_amb(signal, LR, iAngle, level, configurationSetup)
% ICEBERG_SET_AMB  Render signal through Ambisonics decoding with optional
% per-LS calibration applied to the dry signal before convolution.
%
%   signal_ambisonics = iceberg_set_amb(signal, LR, iAngle, level, configurationSetup)

signal_nor = native_normalize_dat(signal);

% Per-LS frequency filter + SPL alignment on the dry signal (matches the
% original ambisonics_set_level positioning: cal then convolve)
if isfield(configurationSetup, 'iLoudspeakerFreqFilter') && ...
   ~isempty(configurationSetup.iLoudspeakerFreqFilter)
    signal_nor = calibrate_ambisonics(signal_nor, level, iAngle, configurationSetup);
end

sinal_Ambisonics_LR = native_convolve(signal_nor, LR);

[D4,~] = ambiDecoder(configurationSetup.ls_dir,'SAD',1,1);

amb_time = decodeBformat(sinal_Ambisonics_LR.time, D4);

signal_ambisonics.time = amb_time;
signal_ambisonics.samplingRate = signal.samplingRate;
signal_ambisonics.nChannels = size(amb_time, 2);
signal_ambisonics.nSamples = size(amb_time, 1);

end
