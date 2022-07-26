function [signalcalibrated] = set_level_ambisonics_fly_in(signal,level,iAngles,configurationSetup)
%%Play and record using the most recent calibration
%Input ItaAudio

new_Level_Factor = configurationSetup.new_Level_Factor;
iFactor = configurationSetup.iFactor;
iLoudspeakerFreqFilter = configurationSetup.iLoudspeakerFreqFilter;
LSArray = configurationSetup.LSArray;

%% Load the most recent calibration
Level_Factor = new_Level_Factor;
iFs = signal.samplingRate;                %Sample Frequency
lsdBperVolt = (20*log10((iFactor)/2e-5));
if level > 95
    error ('Level max is 95 dB');
end


%%
signal_to_play = signal;

%%
activeLSNumbers = zeros(1,length(configurationSetup.ls_dir));
for indexActiveLSNumbers= 1:length(configurationSetup.ls_dir)
    iArrayP = find(LSArray==configurationSetup.ls_dir(indexActiveLSNumbers,1),1);
    activeLSNumbers(indexActiveLSNumbers) = iArrayP;
end
for iCount = activeLSNumbers
    Interpolation(:,iCount) = pchip(iLoudspeakerFreqFilter(iCount).freqVector,...
        iLoudspeakerFreqFilter(iCount).freq,signal_to_play.freqVector);
end

frequencyFilter = itaAudio(Interpolation,iFs,'freq');

% % %%
new_page = zeros(length(signal_to_play.time),1);
stimulus = zeros(length(signal_to_play.time),1);
signal_run = zeros(length(signal_to_play.time),1);



%% Select the filter fo the audio according to the LS (Virtual Loudspeakers are filtered with the Nearest Speaker NSP)
% This is to fit level the ambisonics audio part
if iAngles > 90 && iAngles <= 180
    if  iAngles <= 135
        iChannel = activeLSNumbers(find(configurationSetup.ls_dir==90,1));
    else
        iChannel = activeLSNumbers(find(configurationSetup.ls_dir==180,1));
    end
elseif iAngles > 180 && iAngles <= 270
    if  iAngles <= 225
        iChannel = activeLSNumbers(find(configurationSetup.ls_dir==180,1));
    else
        iChannel = activeLSNumbers(find(configurationSetup.ls_dir==270,1));
    end
elseif iAngles > 270 && iAngles <= 360
    if  iAngles <= 315
        iChannel = activeLSNumbers(find(configurationSetup.ls_dir==270,1));
    else
        iChannel = activeLSNumbers(find(configurationSetup.ls_dir==0,1));
    end
elseif iAngles >= 0 && iAngles <= 90
    if  iAngles <= 45
        iChannel = activeLSNumbers(find(configurationSetup.ls_dir==0,1));
    else
        iChannel = activeLSNumbers(find(configurationSetup.ls_dir==90,1));
    end
end

signal_run_ita_FILTER       = itaAudio(zeros(length(signal_to_play.time),1),iFs,'time');
signal_run_ita_FILTER_LEVEL = itaAudio(zeros(length(signal_to_play.time),1),iFs,'time');

if isnan(signal_to_play.time) 
    new_page(:,1) = 0;
else
    new_page(:,1) = signal_to_play.time;
end
scaler = (sqrt(mean(new_page(:,1).^2)))*2;
stimulus(:,1) = new_page(:,1)./scaler;
if level ~= 'n' %Verify this
    signal_run(:,1) = stimulus(:,1).*repmat(10.^((level - lsdBperVolt )./20),length(new_page(:,1)),1);
else
    signal_run(:,1)= signal_to_play.time;
end
signal_run_ita_dB = itaAudio(signal_run,iFs,'time');
filtered = ita_multiply_spk(signal_run_ita_dB,frequencyFilter.ch(iChannel));
signal_run_ita_FILTER.time(:,1) = filtered.time;
signal_run_ita_FILTER_LEVEL.time(:,1) = signal_run_ita_FILTER.time.*(Level_Factor(iChannel));


signalcalibrated = signal_run_ita_FILTER_LEVEL;


end

