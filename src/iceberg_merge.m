function [final_audio, VBAP_DSER_Part, Amb_LR_Part] = ...
    iceberg_merge(signal, DSER, LR, level, angle, configurationSetup)
activeLSNumbers = zeros(1,length(configurationSetup.ls_dir));

for indexx= 1:length(configurationSetup.ls_dir)
    iArrayP = find(configurationSetup.LSArray==configurationSetup.ls_dir(indexx,1),1);
    activeLSNumbers(indexx) = iArrayP;
end

iFs            = signal.samplingRate;
VBAP_DSER_Part = itaAudio();
Amb_LR_Part  = itaAudio();
VBAP_DSER_Part.samplingRate = iFs;
Amb_LR_Part.samplingRate = iFs;

%% Signal with Omnidirectional reverberation
convolved_DSER_signal = ita_convolve(signal,DSER);
VBAP_DS =   vbap_set_level(convolved_DSER_signal,level,angle,configurationSetup);
VBAP_DS =   VBAP_DS*max(DSER.time); %?

%% fetch VBAP Direct Sound to the Array

for i = 1:length(activeLSNumbers)
    VBAP_DSER_Part.time(:,activeLSNumbers(i)) = VBAP_DS.time(:,i);
end

signal_nor  = ita_normalize_dat(signal,'allchannels','true');
signalAmb   = ambisonics_set_level(signal_nor,level,angle,configurationSetup); 
sinal_Ambisonics_LR   = ita_convolve(signalAmb,LR);         
%% Calibrated Ambisonics to specified Sound Pressure Level
%Decoder Matrix with n=4 LS SAD decoder! %SAD | MMD | EPAD | ALLRAD | CSAD
[D4,~] = ambiDecoder(configurationSetup.ls_dir,'SAD',1,1); 
%Decode to Loudspeaker array
Ambisonics_ERLRSignal = itaAudio(decodeBformat(...
    sinal_Ambisonics_LR.time, D4), iFs, 'time');
%fetch Ambisonics Late reverberation to the Array
for i = 1:length(activeLSNumbers)
    Amb_LR_Part.time(:,activeLSNumbers(i))  = Ambisonics_ERLRSignal.time(:,i);
end
%% Combine DS ER and LR
final_audio = ita_add(VBAP_DSER_Part, Amb_LR_Part);    
if length(configurationSetup.LSArray) > final_audio.dimensions
    final_audio.time(:,final_audio.dimensions+1:length(configurationSetup.LSArray))...
        = zeros(final_audio.nSamples,length(configurationSetup.LSArray)-(final_audio.dimensions));
end

end
