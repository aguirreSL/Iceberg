function [ new_Level_Factor , old_Level_Factor] = getLevel(iFactor,iLoudspeakerFreqFilter,tolerance,increment,nAverage,excitation_signal,iChannel)
[~, playDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('playDeviceID'));
[~, recDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('recDeviceID'));
playDeviceID = playDeviceInfo.deviceID; recDeviceID = recDeviceInfo.deviceID;

upperFreq = 20000;
lowerFreq = 60;

try load('Current_Calibration.mat', 'new_Level_Factor');
catch error (sprintf('Calibration file not found, copy ''current_calibration'' to %s',pwd))
end
old_Level_Factor = new_Level_Factor; %Store last
iFs = 48000;             %Sample Frequency
% tolerance = .5;
% increment = 0.15;
% nAverage = 2;
lsdBperVolt = (20*log10((iFactor)/2e-5));
% max_safe_voltage =  10.^((max_safe_level-lsdBperVolt+20)./20);
if excitation_signal == 1
    signal_to_play = ita_generate('swenlinsweep',[lowerFreq upperFreq],0.1,iFs,16);
elseif excitation_signal == 2
    signal_to_play = ita_generate('pinknoise',1,iFs,16);
    signal_to_play = ita_time_window(signal_to_play,[0.02 0],'time');
    signal_to_play = ita_filter_bandpass(signal_to_play,'upper',upperFreq,'lower',lowerFreq);
    
elseif excitation_signal == 3
    %% LTASS
    % Add your wav file
    file = ita_read(uigetfile({'.\Auxiliar\*.wav'},'Pick a file'));
    % Or you can use this instead:
    % file = ita_read('HINT_Danish_Noise_Female2_Inf_new.wav');
    file = ita_resample(file,iFs);
    userAnswer = questdlg('Would you like crop the silent parts using VAD?', ...
        'Yes, please','No, thank you');
    % VAD please, you should fine tune the VAD accordingly your audio
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
    signal_to_play = ita_time_crop(LTASS,[0 lengthSignal],'time');


end
% Normalize input
signal_to_play = ita_normalize_dat(signal_to_play,'allchannels','true');
for iCount = 1:24
    Interpolation(:,iCount) = pchip(iLoudspeakerFreqFilter(iCount).freqVector,...
        iLoudspeakerFreqFilter(iCount).freq,signal_to_play.freqVector);
end

%Transform filter in itaAudio
% frequencyFilter = ita_time_window((itaAudio(Interpolation,iFs,'freq')),[0 0.4],'time','@hann');
frequencyFilter = (itaAudio(Interpolation,iFs,'freq'));
%%
Level_Factor(1:24) = abs(new_Level_Factor);
% Level_Factor(1:24) = 1;
level = 70;

for iLoudspeaker =1:24
    
%     for   iAverage = 1:nAverage
for iRepeat = 1:nAverage
        new_page = zeros(length(signal_to_play.time),24);
        new_page(:,iLoudspeaker) = signal_to_play.time;
        scaler = (sqrt(mean(new_page(:,iLoudspeaker).^2)))*2;
        stimulus = new_page(:,iLoudspeaker)./scaler; %2.5
        %NPS dB Value Based on dB/V factor
        signal_run = stimulus.*repmat(10.^((level - lsdBperVolt)./20),length(new_page(:,iLoudspeaker)),1);
        %add freq filter to sweep signal
        signal_run_ita_LEVEL = itaAudio(signal_run,iFs,'time');
        selectedFilter = ita_split(frequencyFilter,iLoudspeaker);
        signal_run_ita_FILTER_LEVEL = ita_multiply_spk(signal_run_ita_LEVEL,selectedFilter);
        %Adjust sweep Level
        signal_run_ita_FILTER_LEVEL = signal_run_ita_FILTER_LEVEL.*Level_Factor(iLoudspeaker);
        signal_run_ita_FILTER_LEVEL = ita_time_window(signal_run_ita_FILTER_LEVEL,[0.05 0],'time','windowType', 'hann','symmetric','true');
        signal_run_ita_FILTER_LEVEL = ita_filter_bandpass(signal_run_ita_FILTER_LEVEL,'upper',upperFreq,'lower',lowerFreq);


        new_page(:,iLoudspeaker) = signal_run_ita_FILTER_LEVEL.time;
        %Play 
        playrec('init', iFs, playDeviceID, recDeviceID);
        page = playrec('playrec',new_page,1:24,length(new_page),iChannel);
        while playrec('isFinished')~=1,pause(0.1);end
        recording = playrec('getRec',page);
        recordedSignal = recording;
        
        playrec('reset')
        recorded_signal = itaAudio(recording,iFs,'time');
        recorded_signal = ita_time_window(recorded_signal,[0.05 0],'time','windowType', 'hann','symmetric','true');
        recorded_signal = ita_filter_bandpass(recorded_signal,'upper',upperFreq,'lower',lowerFreq);
        iSignal_Filter_Level = ita_spk2level(recorded_signal*iFactor,0,'added');
        ispl   = 20*log10(iSignal_Filter_Level.freqData/2e-5);
        
        SPLrepeat(iLoudspeaker,iRepeat) = ispl;
        
        SPLaverage(iLoudspeaker) = mean(SPLrepeat(iLoudspeaker,(1:iRepeat)));
        fprintf('  LS = %i n = %i\nSPL %.2f [dB] average SPL %.2f [dB]\n\n',iLoudspeaker,iRepeat,ispl,SPLaverage(iLoudspeaker))
        
        
end        
  
%         while SPLaverage(iLoudspeaker,iAverage)> (level+tolerance)
%             
%             updated_Value = Level_Factor(iLoudspeaker)-increment;
%             
%             Level_Factor(iLoudspeaker) = updated_Value;
%             new_page = zeros(length(signal_to_play.time),24);
%             new_page(:,iLoudspeaker) = signal_to_play.time;
%             scaler = (sqrt(mean(new_page(:,iLoudspeaker).^2)))*2;
%             stimulus = new_page(:,iLoudspeaker)./scaler; %2.5
%             %NPS dB Value Based on dB/V factor
%             signal_run = stimulus.*repmat(10.^((level - lsdBperVolt)./20),length(new_page(:,iLoudspeaker)),1);
%             %add freq filter
%             signal_run_ita_LEVEL = itaAudio(signal_run,iFs,'time');
%             signal_run_ita_FILTER_LEVEL = ita_multiply_spk(signal_run_ita_LEVEL,frequencyFilter.ch(iLoudspeaker));
%             signal_run_ita_FILTER_LEVEL = signal_run_ita_FILTER_LEVEL.*Level_Factor(iLoudspeaker);
%             new_page(:,iLoudspeaker) = signal_run_ita_FILTER_LEVEL.time;
%             playrec('init', iFs, playDeviceID, recDeviceID);
%             page = playrec('playrec',new_page,1:24,length(new_page),iChannel);
%             while playrec('isFinished')~=1,pause(0.1);end
%             recording = playrec('getRec',page);
%             recordedSignal = recording;
%             playrec('reset')
%             recorded_signal = itaAudio(recording,iFs,'time');
%             iSignal_Filter_Level = ita_spk2level(recorded_signal*iFactor,0,'added');
%             ispl   = 20*log10(iSignal_Filter_Level.freqData/2e-5);
%             SPLaverage(iLoudspeaker,iAverage) = ispl;
%             fprintf('  LS = %i\n SPL %.2f dB\n',iLoudspeaker,ispl)
%             %% New up
%             while SPLaverage(iLoudspeaker,iAverage)<(level-tolerance)
%                 Level_Factor(iLoudspeaker) = updated_Value+increment;
%                 new_page = zeros(length(signal_to_play.time),24);
%                 new_page(:,iLoudspeaker) = signal_to_play.time;
%                 scaler = (sqrt(mean(new_page(:,iLoudspeaker).^2)))*2;
%                 stimulus = new_page(:,iLoudspeaker)./scaler; %2.5
%                 %NPS dB Value Based on dB/V factor
%                 signal_run = stimulus.*repmat(10.^((level - lsdBperVolt)./20),length(new_page(:,iLoudspeaker)),1);
%                 %add freq filter
%                 signal_run_ita_LEVEL = itaAudio(signal_run,iFs,'time');
%                 signal_run_ita_FILTER_LEVEL = ita_multiply_spk(signal_run_ita_LEVEL,frequencyFilter.ch(iLoudspeaker));
%                 %Adjust Level
%                 signal_run_ita_FILTER_LEVEL = signal_run_ita_FILTER_LEVEL.*Level_Factor(iLoudspeaker);
%                 new_page(:,iLoudspeaker) = signal_run_ita_FILTER_LEVEL.time;
%                 %%
%                 playrec('init', iFs, playDeviceID, recDeviceID);
%                 page = playrec('playrec',new_page,1:24,length(new_page),iChannel);
%                 while playrec('isFinished')~=1,pause(0.1);end
%                 recording = playrec('getRec',page);
%                 recordedSignal = recording;
%                 playrec('reset')
%                 recorded_signal = itaAudio(recording,iFs,'time');
%                 iSignal_Filter_Level = ita_spk2level(recorded_signal*iFactor,0,'added');
%                 ispl   = 20*log10(iSignal_Filter_Level.freqData/2e-5);
%                 SPLaverage(iLoudspeaker,iAverage) = ispl;
%                 fprintf('  LS = %i\n SPL %.2f dB\n',iLoudspeaker,ispl)
%                 
%             end
%         end
%         while mean(SPLaverage(iLoudspeaker,iAverage))< (level-tolerance)
%             updated_Value = Level_Factor(iLoudspeaker)+increment;
%             Level_Factor(iLoudspeaker) = updated_Value;
%             new_page = zeros(length(signal_to_play.time),24);
%             new_page(:,iLoudspeaker) = signal_to_play.time;
%             scaler = (sqrt(mean(new_page(:,iLoudspeaker).^2)))*2;
%             stimulus = new_page(:,iLoudspeaker)./scaler; %2.5
%             %NPS dB Value Based on dB/V factor
%             signal_run = stimulus.*repmat(10.^((level - lsdBperVolt)./20),length(new_page(:,iLoudspeaker)),1);
%             %add freq filter
%             signal_run_ita_LEVEL = itaAudio(signal_run,iFs,'time');
%             signal_run_ita_FILTER_LEVEL = ita_multiply_spk(signal_run_ita_LEVEL,frequencyFilter.ch(iLoudspeaker));
%             %Adjust Level
%             signal_run_ita_FILTER_LEVEL = signal_run_ita_FILTER_LEVEL.*Level_Factor(iLoudspeaker);
%             new_page(:,iLoudspeaker) = signal_run_ita_FILTER_LEVEL.time;
%             %%
%             playrec('init', iFs, playDeviceID, recDeviceID);
%             page = playrec('playrec',new_page,1:24,length(new_page),iChannel);
%             while playrec('isFinished')~=1,pause(0.1);end
%             recording = playrec('getRec',page);
% %             recordedSignal = recording;
%             playrec('reset')
%             recorded_signal = itaAudio(recording,iFs,'time');
%             iSignal_Filter_Level = ita_spk2level(recorded_signal*iFactor,0,'added');
%             ispl   = 20*log10(iSignal_Filter_Level.freqData/2e-5);
%             SPLaverage(iLoudspeaker,iAverage) = ispl;
%             fprintf('  LS = %i\n SPL %.2f dB\n',iLoudspeaker,ispl)
%             
%         end
%         
%         %% New Down
%         while mean(SPLaverage(iLoudspeaker,iAverage))> (level+tolerance)
%             Level_Factor(iLoudspeaker) = updated_Value-increment;
%             new_page = zeros(length(signal_to_play.time),24);
%             new_page(:,iLoudspeaker) = signal_to_play.time;
%             scaler = (sqrt(mean(new_page(:,iLoudspeaker).^2)))*2;
%             stimulus = new_page(:,iLoudspeaker)./scaler; %2.5
%             %NPS dB Value Based on dB/V factor
%             signal_run = stimulus.*repmat(10.^((level - lsdBperVolt)./20),length(new_page(:,iLoudspeaker)),1);
%             %add freq filter
%             signal_run_ita_LEVEL = itaAudio(signal_run,iFs,'time');
%             signal_run_ita_FILTER_LEVEL = ita_multiply_spk(signal_run_ita_LEVEL,frequencyFilter.ch(iLoudspeaker));
%             signal_run_ita_FILTER_LEVEL = signal_run_ita_FILTER_LEVEL.*Level_Factor(iLoudspeaker);
%             new_page(:,iLoudspeaker) = signal_run_ita_FILTER_LEVEL.time;
%             playrec('init', iFs, playDeviceID, recDeviceID);
%             page = playrec('playrec',new_page,1:24,length(new_page),iChannel);
%             while playrec('isFinished')~=1,pause(0.1);end
%             recording = playrec('getRec',page);
%             recordedSignal = recording;
%             playrec('reset')
%             recorded_signal = itaAudio(recording,iFs,'time');
%             iSignal_Filter_Level = ita_spk2level(recorded_signal*iFactor,0,'added');
%             ispl   = 20*log10(iSignal_Filter_Level.freqData/2e-5);
%             SPLaverage(iLoudspeaker,iAverage) = ispl;
%             fprintf('  LS = %i\n SPL %.2f dB\n',iLoudspeaker,ispl)
%         end
 
while SPLaverage(iLoudspeaker)> (level+tolerance) || SPLaverage(iLoudspeaker)<(level-tolerance)
    
    if SPLaverage(iLoudspeaker)>(level+tolerance)
        for iRepeat = 1:nAverage
        updated_Value = Level_Factor(iLoudspeaker)-increment;
        Level_Factor(iLoudspeaker) = updated_Value;
        new_page = zeros(length(signal_to_play.time),24);
        new_page(:,iLoudspeaker) = signal_to_play.time;
        scaler = (sqrt(mean(new_page(:,iLoudspeaker).^2)))*2;
        stimulus = new_page(:,iLoudspeaker)./scaler; %2.5
        %NPS dB Value Based on dB/V factor
        signal_run = stimulus.*repmat(10.^((level - lsdBperVolt)./20),length(new_page(:,iLoudspeaker)),1);
        %add freq filter
        signal_run_ita_LEVEL = itaAudio(signal_run,iFs,'time');
        selectedFilter = ita_split(frequencyFilter,iLoudspeaker);
        signal_run_ita_FILTER_LEVEL = ita_multiply_spk(signal_run_ita_LEVEL,selectedFilter);
        signal_run_ita_FILTER_LEVEL = signal_run_ita_FILTER_LEVEL.*Level_Factor(iLoudspeaker);
        signal_run_ita_FILTER_LEVEL = ita_filter_bandpass(signal_run_ita_FILTER_LEVEL,'upper',upperFreq,'lower',lowerFreq);

        new_page(:,iLoudspeaker) = signal_run_ita_FILTER_LEVEL.time;
        playrec('init', iFs, playDeviceID, recDeviceID);
        page = playrec('playrec',new_page,1:24,length(new_page),iChannel);
        while playrec('isFinished')~=1,pause(0.1);end
        recording = playrec('getRec',page);
        recordedSignal = recording;
        playrec('reset')
        recorded_signal = itaAudio(recording,iFs,'time');
        recorded_signal = ita_filter_bandpass(recorded_signal,'upper',upperFreq,'lower',lowerFreq);
        iSignal_Filter_Level = ita_spk2level(recorded_signal*iFactor,0,'added');
        ispl   = 20*log10(iSignal_Filter_Level.freqData/2e-5);
        
        SPLrepeat(iLoudspeaker,iRepeat) = ispl;
        SPLaverage(iLoudspeaker) = mean(SPLrepeat(iLoudspeaker,(1:iRepeat)));
        fprintf('  LS = %i n = %i\nSPL %.2f [dB] average SPL %.2f [dB]\n\n',iLoudspeaker,iRepeat,ispl,SPLaverage(iLoudspeaker))
        end  
    elseif SPLaverage(iLoudspeaker)<(level-tolerance)
        for iRepeat = 1:nAverage
        updated_Value = Level_Factor(iLoudspeaker)+increment;
        Level_Factor(iLoudspeaker) = updated_Value;
        new_page = zeros(length(signal_to_play.time),24);
        new_page(:,iLoudspeaker) = signal_to_play.time;
        scaler = (sqrt(mean(new_page(:,iLoudspeaker).^2)))*2;
        stimulus = new_page(:,iLoudspeaker)./scaler; %2.5
        %NPS dB Value Based on dB/V factor
        signal_run = stimulus.*repmat(10.^((level - lsdBperVolt)./20),length(new_page(:,iLoudspeaker)),1);
        %add freq filter
        signal_run_ita_LEVEL = itaAudio(signal_run,iFs,'time');
        selectedFilter = ita_split(frequencyFilter,iLoudspeaker);
        signal_run_ita_FILTER_LEVEL = ita_multiply_spk(signal_run_ita_LEVEL,selectedFilter);
        signal_run_ita_FILTER_LEVEL = signal_run_ita_FILTER_LEVEL.*Level_Factor(iLoudspeaker);
        signal_run_ita_FILTER_LEVEL = ita_filter_bandpass(signal_run_ita_FILTER_LEVEL,'upper',upperFreq,'lower',lowerFreq);

        new_page(:,iLoudspeaker) = signal_run_ita_FILTER_LEVEL.time;
        playrec('init', iFs, playDeviceID, recDeviceID);
        page = playrec('playrec',new_page,1:24,length(new_page),iChannel);
        while playrec('isFinished')~=1,pause(0.1);end
        recording = playrec('getRec',page);
        recordedSignal = recording;
        playrec('reset')
        recorded_signal = itaAudio(recording,iFs,'time');
        recorded_signal = ita_filter_bandpass(recorded_signal,'upper',upperFreq,'lower',lowerFreq);
        iSignal_Filter_Level = ita_spk2level(recorded_signal*iFactor,0,'added');
        ispl   = 20*log10(iSignal_Filter_Level.freqData/2e-5);
        SPLrepeat(iLoudspeaker,iRepeat) = ispl;
        SPLaverage(iLoudspeaker) = mean(SPLrepeat(iLoudspeaker,(1:iRepeat)));
        fprintf('  LS = %i n = %i\nSPL %.2f [dB] average SPL %.2f [dB]\n\n',iLoudspeaker,iRepeat,ispl,SPLaverage(iLoudspeaker))
        end  
    end
end

pause(.5)
        
        
    end
    
    
    
% end
new_Level_Factor = Level_Factor;
end
