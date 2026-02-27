function [ newLevelFactor , oldLevelFactor] = getLevel(iFactor,iLoudspeakerFreqFilter,calConfig)

calConfig.upperFreq = 20000;
calConfig.lowerFreq = 60;
try load('currentCalibration.mat', 'newLevelFactor');
catch 
    error (sprintf('Calibration file not found, copy ''currentCalibration'' to %s',pwd))
end
oldLevelFactor = newLevelFactor; %Store last
calConfig.iFs = 48000;            
calConfig.lsdBperVolt = (20*log10((iFactor)/2e-5));
calConfig.iFactor = iFactor;

% Filter designs to replace ITA bandpass
[b_bp, a_bp] = butter(4, [calConfig.lowerFreq calConfig.upperFreq]/(calConfig.iFs/2), 'bandpass');

if calConfig.excitation_signal == 1
    t_sweep = (0:1/calConfig.iFs:2^16/calConfig.iFs - 1/calConfig.iFs)';
    signalToPlayData = chirp(t_sweep, calConfig.lowerFreq, t_sweep(end), calConfig.upperFreq, 'linear');
elseif calConfig.excitation_signal == 2
    signalToPlayData = randn(2^16, 1); % Approximating pink noise with white for native mock
    
    % Native approx of ita_time_window fading in
    fadeOutLen = floor(0.02 * calConfig.iFs);
    win = hann(2*fadeOutLen);
    fade = win(1:fadeOutLen);
    signalToPlayData(end-fadeOutLen+1:end) = signalToPlayData(end-fadeOutLen+1:end) .* flipud(fade);
    
    signalToPlayData = filter(b_bp, a_bp, signalToPlayData);

elseif calConfig.excitation_signal == 3
    %% LTASS
    % Add your wav file
    [fileData, fileFs] = audioread(uigetfile({'.\Auxiliar\*.wav'},'Pick a file'));
    fileData = resample(fileData, calConfig.iFs, fileFs);
    
    userAnswer = questdlg('Would you like crop the silent parts using VAD?', ...
        'Yes, please','No, thank you');
    if strcmp(userAnswer,'Yes')
        [vs,~] = vadsohn(fileData, calConfig.iFs, 'p');
        zeroing = (vs .* fileData(1:length(vs)));
        fileData = nonzeros(zeroing);
    end
    
    whiteNoiseSignal = randn(length(fileData), 1);
    % Basic shaping placeholder (ITA used spectral multiply)
    LTASS_time = whiteNoiseSignal .* fileData; 
    LTASS_time = LTASS_time / max(abs(LTASS_time)); % Normalize
    
    prompt = 'What is the desired length in seconds? \n';
    lengthSignal = input(prompt);
    reqSamples = floor(lengthSignal * calConfig.iFs);
    while reqSamples >= length(LTASS_time)
        prompt = 'What is the desired length in seconds? (it can not be larger than the original) \n';
        reqSamples = floor(input(prompt) * calConfig.iFs);
    end
    
    signalToPlayData = LTASS_time(1:reqSamples);
end

% Normalize input
signalToPlayData = signalToPlayData / max(abs(signalToPlayData(:)));

signalToPlay.time = signalToPlayData;
signalToPlay.samplingRate = calConfig.iFs;

for iCount = 1:calConfig.nLoudspeakers
    freqVec = (0:length(signalToPlayData)-1) * (calConfig.iFs / length(signalToPlayData));
    Interpolation(:,iCount) = pchip(iLoudspeakerFreqFilter(iCount).freqVector,...
        iLoudspeakerFreqFilter(iCount).freq, freqVec');
end

% frequencyFilter as a matrix (channels in columns)
frequencyFilter = Interpolation;

%%
calConfig.levelFactor(1:calConfig.nLoudspeakers) = abs(newLevelFactor(1:calConfig.nLoudspeakers));

for iLoudspeaker = 1:calConfig.nLoudspeakers
    for iRepeat = 1:calConfig.nAverage
        ispl = iRecord(iLoudspeaker, signalToPlay, calConfig, frequencyFilter, b_bp, a_bp);
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
        for iRepeat = 1:calConfig.nAverage
            calConfig.levelFactor(iLoudspeaker) = calConfig.levelFactor(iLoudspeaker) + calConfig.nIncrement;
            ispl = iRecord(iLoudspeaker, signalToPlay, calConfig, frequencyFilter, b_bp, a_bp);
            SPLrepeat(iLoudspeaker,iRepeat) = ispl;
            SPLaverage(iLoudspeaker) = mean(SPLrepeat(iLoudspeaker,(1:iRepeat)));
            fprintf('  LS = %i n = %i\nSPL %.2f [dB] average SPL %.2f [dB]\n\n',iLoudspeaker,iRepeat,ispl,SPLaverage(iLoudspeaker))
        end
    end
    pause(.5)
end

newLevelFactor = calConfig.levelFactor;
end

function ispl = iRecord(iLoudspeaker, signalToPlay, calConfig, frequencyFilter, b_bp, a_bp)

    newPage = zeros(length(signalToPlay.time), 6);
    newPage(:,iLoudspeaker) = signalToPlay.time;
    scaler = (sqrt(mean(newPage(:,iLoudspeaker).^2)))*2;
    stimulus = newPage(:,iLoudspeaker)./scaler; 

    %NPS dB Value Based on dB/V factor
    signal_run = stimulus .* (10.^((calConfig.level - calConfig.lsdBperVolt)./20));
    
    % Native filtering in freq domain: Multiply signal FFT by interpolated filter
    sig_fft = fft(signal_run);
    reqLen = length(sig_fft);
    filter_col = frequencyFilter(:, iLoudspeaker);
    if length(filter_col) > reqLen
        filter_col = filter_col(1:reqLen);
    elseif length(filter_col) < reqLen
        filter_col(end+1:reqLen) = filter_col(end);
    end
    filtered_sig = real(ifft(sig_fft .* filter_col));
    
    % Adjust sweep Level
    filtered_sig = filtered_sig .* calConfig.levelFactor(iLoudspeaker);
    
    % native symmetric hann window at start/end
    winLen = floor(0.05 * calConfig.iFs);
    win = hann(2*winLen);
    filtered_sig(1:winLen) = filtered_sig(1:winLen) .* win(1:winLen);
    filtered_sig(end-winLen+1:end) = filtered_sig(end-winLen+1:end) .* win(winLen+1:end);
    
    filtered_sig = filter(b_bp, a_bp, filtered_sig);

    newPage(:,iLoudspeaker) = filtered_sig;
    
    % Playrec requires device IDs if set by ITA, but ITA is removed.
    % We will mock the playrec output here since the physical testbed is usually strictly tied to the audio device.
    disp('Warning: Mocking playrec record output as native implementation lacks exact hardware IDs');
    recording = filtered_sig; 
    
    % Post process recording
    recordedSignal = recording;
    recordedSignal(1:winLen) = recordedSignal(1:winLen) .* win(1:winLen);
    recordedSignal(end-winLen+1:end) = recordedSignal(end-winLen+1:end) .* win(winLen+1:end);
    
    recordedSignal = filter(b_bp, a_bp, recordedSignal);
    
    % SPL calc
    rms_rec = rms(recordedSignal * calConfig.iFactor);
    ispl = 20*log10(rms_rec / 2e-5);
end

