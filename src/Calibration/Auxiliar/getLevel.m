function [ newLevelFactor , oldLevelFactor] = getLevel(iFactor,iLoudspeakerFreqFilter,calConfig)
[~, playDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('playDeviceID'));
[~, recDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('recDeviceID'));
calConfig.playDeviceID = playDeviceInfo.deviceID;
calConfig.recDeviceID = recDeviceInfo.deviceID;

calConfig.upperFreq = 20000;
calConfig.lowerFreq = 60;
try load('currentCalibration.mat', 'newLevelFactor');
catch error (sprintf('Calibration file not found, copy ''currentCalibration'' to %s',pwd))
end
oldLevelFactor = newLevelFactor; %Store last
calConfig.iFs = 48000;            
calConfig.lsdBperVolt = (20*log10((iFactor)/2e-5));
calConfig.iFactor = iFactor;
% max_safe_voltage =  10.^((max_safe_level-lsdBperVolt+20)./20);

if calConfig.excitation_signal == 1
    signalToPlay = ita_generate('swenlinsweep',[calConfig.lowerFreq calConfig.upperFreq],0.1,calConfig.iFs,16);
elseif calConfig.excitation_signal == 2
    signalToPlay = ita_generate('pinknoise',1,calConfig.iFs,16);
    signalToPlay = ita_time_window(signalToPlay,[0.02 0],'time');
    signalToPlay = ita_filter_bandpass(signalToPlay,'upper',calConfig.upperFreq,'lower',calConfig.lowerFreq);

elseif calConfig.excitation_signal == 3
    %% LTASS
    % Add your wav file
    file = ita_read(uigetfile({'.\Auxiliar\*.wav'},'Pick a file'));
    % Or you can use this instead:
    % file = ita_read('HINT_Danish_Noise_Female2_Inf_new.wav');
    file = ita_resample(file,iFs);
    userAnswer = questdlg('Would you like crop the silent parts using VAD?', ...
        'Yes, please','No, thank you');
    % VAD please, you can fine tune the VAD to your audio
    if strcmp(userAnswer,'Yes')
        [vs,~]=vadsohn(file.time,iFs,'p');
        zeroing = (vs'.*file.time(1:length(vs))');
        file = itaAudio(nonzeros(zeroing),iFs,'time');
    end
    smoothSignal        = ita_smooth_frequency(file);
    whiteNoiseSignal    = ita_generate('noise',1,file.samplingRate,file.fftDegree);
    % pink            = ita_generate('pinknoise',1,osinal.samplingRate,lengthSignal);
    LTASS            = ita_multiply_spk(whiteNoiseSignal,smoothSignal);
    LTASS            = ita_normalize_dat(LTASS);
    %%
    prompt = 'What is the desired length in seconds? \n';
    lengthSignal = input(prompt);
    while lengthSignal >= LTASS.trackLength
        prompt = 'What is the desired length in seconds? (it can not be larger than the original) \n';
        lengthSignal = input(prompt);
    end
    signalToPlay = ita_time_crop(LTASS,[0 lengthSignal],'time');


end
% Normalize input
signalToPlay = ita_normalize_dat(signalToPlay,'allchannels','true');
for iCount = 1:calConfig.nLoudspeakers
    Interpolation(:,iCount) = pchip(iLoudspeakerFreqFilter(iCount).freqVector,...
        iLoudspeakerFreqFilter(iCount).freq,signalToPlay.freqVector);
end

%Transform filter in itaAudio
% frequencyFilter = ita_time_window((itaAudio(Interpolation,iFs,'freq')),[0 0.4],'time','@hann');
frequencyFilter = (itaAudio(Interpolation,calConfig.iFs,'freq'));
%%
calConfig.levelFactor(1:calConfig.nLoudpeakers) = abs(newLevelFactor(1:calConfig.nLoudpeakers));

for iLoudspeaker = 1:calConfig.nLoudpeakers

    for iRepeat = 1:calConfig.nAverage

        ispl = iRecord(iLoudspeaker,signalToPlay,calConfig,frequencyFilter);
        SPLrepeat(iLoudspeaker,iRepeat) = ispl;
        SPLaverage(iLoudspeaker) = mean(SPLrepeat(iLoudspeaker,(1:iRepeat)));
        fprintf('  LS = %i n = %i\nSPL %.2f [dB] average SPL %.2f [dB]\n\n',iLoudspeaker,iRepeat,ispl,SPLaverage(iLoudspeaker))

    end


    while SPLaverage(iLoudspeaker) > (calConfig.level + calConfig.nTolerance) || SPLaverage(iLoudspeaker) < (calConfig.level - calConfig.nTolerance)
        if SPLaverage(iLoudspeaker) > (calConfig.level + calConfig.nTolerance)
            calConfig.nIncrement = -abs(calConfig.nIncrement); % Decrease level
        else
            calConfig.nIncrement = abs(calConfig.nIncrement); % Increase level
        end
        %         calConfig.levelFactor(iLoudspeaker) = calConfig.levelFactor(iLoudspeaker) + calConfig.nIncrement;
        for iRepeat = 1:calConfig.nAverage
            calConfig.levelFactor(iLoudspeaker) = calConfig.levelFactor(iLoudspeaker) + calConfig.nIncrement;
            ispl = iRecord(iLoudspeaker,signalToPlay,calConfig,frequencyFilter);
            SPLrepeat(iLoudspeaker,iRepeat) = ispl;
            SPLaverage(iLoudspeaker) = mean(SPLrepeat(iLoudspeaker,(1:iRepeat)));
            fprintf('  LS = %i n = %i\nSPL %.2f [dB] average SPL %.2f [dB]\n\n',iLoudspeaker,iRepeat,ispl,SPLaverage(iLoudspeaker))
        end
    end

    pause(.5)

end

newLevelFactor = calConfig.levelFactor;
end


function ispl = iRecord(iLoudspeaker,signalToPlay,calConfig,frequencyFilter)

newPage = zeros(length(signalToPlay.time),6);
newPage(:,iLoudspeaker) = signalToPlay.time;
scaler = (sqrt(mean(newPage(:,iLoudspeaker).^2)))*2;
stimulus = newPage(:,iLoudspeaker)./scaler; %2.5
%NPS dB Value Based on dB/V factor
signal_run = stimulus.*repmat(10.^((calConfig.level - calConfig.lsdBperVolt)./20),length(newPage(:,iLoudspeaker)),1);
%add freq filter to sweep signal
signal_run_ita_LEVEL = itaAudio(signal_run,calConfig.iFs,'time');
selectedFilter = ita_split(frequencyFilter,iLoudspeaker);
signal_run_ita_FILTER_LEVEL = ita_multiply_spk(signal_run_ita_LEVEL,selectedFilter);
%Adjust sweep Level
signal_run_ita_FILTER_LEVEL = signal_run_ita_FILTER_LEVEL.*calConfig.levelFactor(iLoudspeaker);
signal_run_ita_FILTER_LEVEL = ita_time_window(signal_run_ita_FILTER_LEVEL,[0.05 0],'time','windowType', 'hann','symmetric','true');
signal_run_ita_FILTER_LEVEL = ita_filter_bandpass(signal_run_ita_FILTER_LEVEL,'upper',calConfig.upperFreq,'lower',calConfig.lowerFreq);

newPage(:,iLoudspeaker) = signal_run_ita_FILTER_LEVEL.time;
%Play
playrec('init', calConfig.iFs, calConfig.playDeviceID, calConfig.recDeviceID);
page = playrec('playrec',newPage,1:calConfig.nLoudpeakers,length(newPage),calConfig.iChannel);
while playrec('isFinished')~=1,pause(0.1);end
recording = playrec('getRec',page);
%         recordedSignal = recording;

playrec('reset')
recordedSignal = itaAudio(recording,calConfig.iFs,'time');
recordedSignal = ita_time_window(recordedSignal,[0.05 0],'time','windowType', 'hann','symmetric','true');
recordedSignal = ita_filter_bandpass(recordedSignal,'upper',calConfig.upperFreq,'lower',calConfig.lowerFreq);
iSignalFilterLevel = ita_spk2level(recordedSignal*calConfig.iFactor,0,'added');
ispl   = 20*log10(iSignalFilterLevel.freqData/2e-5);

end

