function [final_audio, VBAP_DSER_Part, Amb_LR_Part] = iceberg_merge(signal, DSER, LR, level, angle, configurationSetup)
% Native Implementation
activeLSNumbers = zeros(1,length(configurationSetup.ls_dir));

for indexx = 1:length(configurationSetup.ls_dir)
    iArrayP = find(configurationSetup.lsArray==configurationSetup.ls_dir(indexx,1),1);
    activeLSNumbers(indexx) = iArrayP;
end

iFs = signal.samplingRate;

VBAP_DSER_Part.samplingRate = iFs;
Amb_LR_Part.samplingRate = iFs;

%% Signal with Omnidirectional reverberation
VBAP_DS = iceberg_set_vbap(signal, DSER, angle, configurationSetup);

% Scale by max DSER time peak natively
VBAP_DS.time = VBAP_DS.time * max(max(abs(DSER.time))); 

%% fetch VBAP Direct Sound to the Array
VBAP_DSER_Part.time = zeros(size(VBAP_DS.time, 1), max(activeLSNumbers));
for i = 1:length(activeLSNumbers)
    VBAP_DSER_Part.time(:,activeLSNumbers(i)) = VBAP_DS.time(:,i);
end
VBAP_DSER_Part.nChannels = size(VBAP_DSER_Part.time, 2);
VBAP_DSER_Part.nSamples = size(VBAP_DSER_Part.time, 1);

signalAmb = iceberg_set_amb(signal, LR, configurationSetup);         

%% Calibrated Ambisonics to specified Sound Pressure Level
%Decoder Matrix with n=4 LS SAD decoder! %SAD | MMD | EPAD | ALLRAD | CSAD
[D4,~] = ambiDecoder(configurationSetup.ls_dir,'SAD',1,1); 

%Decode to Loudspeaker array natively
ambDecodedTime = decodeBformat(signalAmb.time, D4);

%fetch Ambisonics Late reverberation to the Array
Amb_LR_Part.time = zeros(size(ambDecodedTime, 1), max(activeLSNumbers));
for i = 1:length(activeLSNumbers)
    Amb_LR_Part.time(:,activeLSNumbers(i))  = ambDecodedTime(:,i);
end
Amb_LR_Part.nChannels = size(Amb_LR_Part.time, 2);
Amb_LR_Part.nSamples = size(Amb_LR_Part.time, 1);

%% Combine DS ER and LR
final_audio = native_add(VBAP_DSER_Part, Amb_LR_Part);    
if length(configurationSetup.lsArray) > final_audio.nChannels
    final_audio.time(:,final_audio.nChannels+1:length(configurationSetup.lsArray))...
        = zeros(final_audio.nSamples,length(configurationSetup.lsArray)-(final_audio.nChannels));
    final_audio.nChannels = size(final_audio.time, 2);
end

end
