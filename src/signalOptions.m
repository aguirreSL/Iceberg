function [selectedSignal, VBAP_placeholder, Amb_placeholder] = signalOptions(Play_This_Signal)
%% Signals (Here one can find some examples, you can define/import any signal) %%% 
sampleFrequency = 44100; 
lengthSignalSeconds = 3; % Seconds
fftDegree = log2(lengthSignalSeconds*sampleFrequency); %From seconds to fft degree
frequencyLimits = [20 20000]; %Check Frequency response from your loudspeakers

switch Play_This_Signal
    case 1
        % signal_1 White Noise
        selectedSignal = ita_generate('whitenoise',1,sampleFrequency,fftDegree);        
    case 2
        % signal_2 Pink Noise
        selectedSignal = ita_generate('pinknoise',1,sampleFrequency,fftDegree);
    case 3
        %signal_3 LTASS - long-term average speech spectrum
        %     file = ita_read('HINT_Danish_Noise_Male1_Inf_new.wav');                       %importing the MALE speech signal 
            file = ita_read('.\wavFiles\HINT_Danish_Noise_Female1_Inf_new.wav');             %importing the FEMALE speech signal
            [vs,~]=vadsohn(file.time,file.samplingRate,'p');                                %Voice Active Detector
            zeroing = [(vs'.*file.time(1:length(vs))'),0];                                  %Mute silent parts  
            file = itaAudio([nonzeros(zeroing);0],file.samplingRate,'time');                %ITA Audio
            smoothSignal        = ita_smooth_frequency(file);                             %Smoth signal (Frequency shape)   
            whiteNoiseSignal    = ita_generate('noise',1,file.samplingRate,file.fftDegree); %Generate white-noise signal/ (is pink more appropriated?)
            LTASS            = ita_multiply_spk(whiteNoiseSignal,smoothSignal);             %Shaping random signal to speech
            LTASS            = ita_time_crop(LTASS,[0 lengthSignalSeconds],'time');         %Defined time (it should be smaller than the imported speech audio)
            LTASS            = ita_normalize_dat(LTASS);                                    %Normalize time domain
            selectedSignal = ita_resample(LTASS,sampleFrequency);                                 %Resampling if needed
    case 4
        % signal_4 Log sweep
        selectedSignal = ita_generate('ccxsweep',frequencyLimits,sampleFrequency,fftDegree);
    case 5
        % signal_5 Linear sweep
        selectedSignal = ita_generate('swenlinsweep',frequencyLimits,0.0,sampleFrequency,fftDegree);
    case 6
        % signal_6 ISTS International Speech Test Signal
        selectedSignal = ita_read([pwd '\wavFiles\ISTS-V1.0_60s_16bit.wav']);
    case 7
        %Pure tone 400 Hz 
        selectedSignal = ita_generate('sine',1,400,sampleFrequency,fftDegree);
end
VBAP_placeholder = itaAudio();
VBAP_placeholder.samplingRate = sampleFrequency;
Amb_placeholder  = itaAudio();
Amb_placeholder.samplingRate = sampleFrequency;
end

