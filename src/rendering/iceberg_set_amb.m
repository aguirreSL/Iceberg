function signal_ambisonics = iceberg_set_amb(signal, LR, configurationSetup)
% Updated: Native Replacement

% Native signal normalization
max_val = max(abs(signal.time(:)));
if max_val > 0
    signal_nor.time = signal.time / max_val;
else
    signal_nor.time = signal.time;
end
signal_nor.samplingRate = signal.samplingRate;
signal_nor.nChannels = size(signal_nor.time, 2);
signal_nor.nSamples = size(signal_nor.time, 1);

% Native signal convolution
sinal_Ambisonics_LR = native_convolve(signal_nor, LR);

% Decoder Matrix with n=4 LS SAD decoder! %SAD | MMD | EPAD | ALLRAD | CSAD
[D4,~] = ambiDecoder(configurationSetup.ls_dir,'SAD',1,1);

% Decode to Loudspeaker array
amb_time = decodeBformat(sinal_Ambisonics_LR.time, D4);

signal_ambisonics.time = amb_time;
signal_ambisonics.samplingRate = signal.samplingRate;
signal_ambisonics.nChannels = size(amb_time, 2);
signal_ambisonics.nSamples = size(amb_time, 1);

end
