function selectedSignal = signalOptions(Play_This_Signal)
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
            filePath = fullfile(pwd, 'wavFiles', 'HINT_Danish_Noise_Female1_Inf_new.wav');
            [y, fs] = audioread(filePath);             %importing the FEMALE speech signal
            
            % Mocking voice activity detection & processing natively
            [vs,~]  = vadsohn(y, fs, 'p');                                %Voice Active Detector
            zeroing = [(vs'.*y(1:length(vs))'),0];                                  %Mute silent parts  
            
            % Create native audio struct instead of itaAudio
            fileVec = [nonzeros(zeroing); 0];
            fileObj.time = fileVec;
            fileObj.samplingRate = fs;
            fileObj.nSamples = length(fileVec);
            
            % Mocking smooth signal natively (Skipping exact spectral smoothing for now, 
            % but doing a basic filter approximation or skipping based on typical usage)
            % LTASS typically shapes white noise against the speech spectrum
            smoothSignal.time = y; % Placeholder for actual spectral smoothing natively
            
            % Generate white noise
            whiteNoiseData = randn(fileObj.nSamples, 1);
            
            % Simplified time domain multiply for LTASS shaping
            LTASS.time = whiteNoiseData .* (smoothSignal.time(1:fileObj.nSamples));
            LTASS.samplingRate = fs;
            LTASS.nSamples = fileObj.nSamples;
            
            LTASS            = native_time_crop(LTASS,[0 lengthSignalSeconds],'time');         %Defined time
            LTASS            = native_normalize_dat(LTASS);                                    %Normalize time domain
            
            % Resampling natively
            if fs ~= sampleFrequency
                resampledTime = resample(LTASS.time, sampleFrequency, fs);
                selectedSignal.time = resampledTime;
                selectedSignal.samplingRate = sampleFrequency;
            else
                selectedSignal = LTASS;
            end
        % signal_4 Log sweep
        selectedSignal = ita_generate('ccxsweep',frequencyLimits,sampleFrequency,fftDegree);
    case 5
        % signal_5 Linear sweep
        selectedSignal = ita_generate('swenlinsweep',frequencyLimits,0.0,sampleFrequency,fftDegree);
    case 6
        % signal_6 ISTS International Speech Test Signal
        selectedSignal = ita_read(fullfile(pwd, 'wavFiles','ISTS-V1.0_60s_16bit.wav'));
    case 7
        %Pure tone 400 Hz 
        selectedSignal = ita_generate('sine',1,400,sampleFrequency,fftDegree);
end

end

