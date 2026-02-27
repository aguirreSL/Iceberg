function selectedSignal = signalOptions(Play_This_Signal)
%% Signals (Here one can find some examples, you can define/import any signal) %%% 
sampleFrequency = 44100; 
lengthSignalSeconds = 3; % Seconds
fftDegree = log2(lengthSignalSeconds*sampleFrequency); %From seconds to fft degree
frequencyLimits = [20 20000]; %Check Frequency response from your loudspeakers

switch Play_This_Signal
    case 1
        % signal_1 White Noise
        nSamples = round(lengthSignalSeconds * sampleFrequency);
        selectedSignal.time = randn(nSamples, 1);
        selectedSignal.samplingRate = sampleFrequency;
        selectedSignal.nSamples = nSamples;
        selectedSignal = native_normalize_dat(selectedSignal);
    case 2
        % signal_2 Pink Noise (1/f spectral shaping of white noise)
        nSamples = round(lengthSignalSeconds * sampleFrequency);
        whiteNoise = randn(nSamples, 1);
        % Apply 1/f filter in frequency domain
        N = length(whiteNoise);
        freqBins = (1:floor(N/2))';
        pinkFilter = 1 ./ sqrt(freqBins);
        X = fft(whiteNoise);
        X(2:floor(N/2)+1) = X(2:floor(N/2)+1) .* [pinkFilter; pinkFilter(end)];
        X(floor(N/2)+2:end) = conj(flipud(X(2:floor(N/2))));
        pinkNoise = real(ifft(X));
        selectedSignal.time = pinkNoise;
        selectedSignal.samplingRate = sampleFrequency;
        selectedSignal.nSamples = nSamples;
        selectedSignal = native_normalize_dat(selectedSignal);
    case 3
        %signal_3 LTASS - long-term average speech spectrum
            projectRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            filePath = fullfile(projectRoot, 'src', 'wavFiles', 'HINT_Danish_Noise_Female1_Inf_new.wav');
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
    case 4
        % signal_4 Log sweep (native chirp replacement for ita_generate ccxsweep)
        nSamples = round(lengthSignalSeconds * sampleFrequency);
        t = (0:nSamples-1)' / sampleFrequency;
        selectedSignal.time = chirp(t, frequencyLimits(1), lengthSignalSeconds, frequencyLimits(2), 'logarithmic');
        selectedSignal.samplingRate = sampleFrequency;
        selectedSignal.nSamples = nSamples;
    case 5
        % signal_5 Linear sweep (native chirp replacement for ita_generate swenlinsweep)
        nSamples = round(lengthSignalSeconds * sampleFrequency);
        t = (0:nSamples-1)' / sampleFrequency;
        selectedSignal.time = chirp(t, frequencyLimits(1), lengthSignalSeconds, frequencyLimits(2), 'linear');
        selectedSignal.samplingRate = sampleFrequency;
        selectedSignal.nSamples = nSamples;
    case 6
        % signal_6 ISTS International Speech Test Signal (native audioread)
        projectRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
        istPath = fullfile(projectRoot, 'src', 'wavFiles','ISTS-V1.0_60s_16bit.wav');
        [istData, istFs] = audioread(istPath);
        selectedSignal.time = istData;
        selectedSignal.samplingRate = istFs;
        selectedSignal.nSamples = size(istData, 1);
    case 7
        % Pure tone 400 Hz
        nSamples = round(lengthSignalSeconds * sampleFrequency);
        t = (0:nSamples-1)' / sampleFrequency;
        selectedSignal.time = sin(2 * pi * 400 * t);
        selectedSignal.samplingRate = sampleFrequency;
        selectedSignal.nSamples = nSamples;
end

end

