function signal_ambisonics = iceberg_set_amb(signal, LR, configurationSetup)
% Updated: 09/06/2025
signal_nor  = ita_normalize_dat(signal,'allchannels','true');
% signalAmb   = ambisonics_set_level(signal_nor,level,angle,configurationSetup);
sinal_Ambisonics_LR   = ita_convolve(signal_nor,LR);
% Decoder Matrix with n=4 LS SAD decoder! %SAD | MMD | EPAD | ALLRAD | CSAD
[D4,~] = ambiDecoder(configurationSetup.ls_dir,'SAD',1,1);
% Decode to Loudspeaker array
signal_ambisonics = itaAudio(decodeBformat(...
    sinal_Ambisonics_LR.time, D4), signal.samplingRate, 'time');

end
