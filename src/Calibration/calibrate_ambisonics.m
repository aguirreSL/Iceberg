function signalcalibrated = calibrate_ambisonics(signal,level,iAngles,configurationSetup)

iLoudspeakerFreqFilter = configurationSetup.iLoudspeakerFreqFilter;
signal_to_play = signal;
new_page = zeros(length(signal_to_play.time),1);
stimulus = zeros(length(signal_to_play.time),1);
signal_run = zeros(length(signal_to_play.time),1);
iFs = signal.samplingRate;                

lsdBperVolt = (20*log10((configurationSetup.iFactor)/2e-5));
if level > 95
    error ('Level max is 95 dB');
end

activeLSNumbers = zeros(1,length(configurationSetup.ls_dir));

for indexActiveLSNumbers= 1:length(configurationSetup.ls_dir)
    iArrayP = find(configurationSetup.LSArray==configurationSetup.ls_dir(indexActiveLSNumbers,1),1);
    activeLSNumbers(indexActiveLSNumbers) = iArrayP;
end

for iCount = activeLSNumbers
    Interpolation(:,iCount) = pchip(iLoudspeakerFreqFilter(iCount).freqVector,...
        iLoudspeakerFreqFilter(iCount).freq,signal_to_play.freqVector);
end

% Native struct for frequency filter (replaces itaAudio freq-domain constructor)
frequencyFilter.freq = Interpolation;
frequencyFilter.samplingRate = iFs;

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

% Native struct to hold filtered/leveled output (replaces itaAudio time-domain constructor)
signal_run_FILTER.time = zeros(length(signal_to_play.time), 1);
signal_run_FILTER.samplingRate = iFs;
signal_run_FILTER_LEVEL.time = zeros(length(signal_to_play.time), 1);
signal_run_FILTER_LEVEL.samplingRate = iFs;

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
% Native spectral multiplication (replaces ita_multiply_spk)
% Multiply signal spectrum by frequency filter in freq domain
sigRunFFT = fft(signal_run);
filterResp = frequencyFilter.freq(:, iChannel);
if length(filterResp) ~= length(sigRunFFT)
    filterResp = interp1(linspace(0,1,length(filterResp)), filterResp, linspace(0,1,length(sigRunFFT))).';
end
filtered = real(ifft(sigRunFFT .* filterResp));
signal_run_FILTER.time(:,1) = filtered;
signal_run_FILTER_LEVEL.time(:,1) = signal_run_FILTER.time .* (configurationSetup.newLevelFactor(iChannel));


signalcalibrated = signal_run_FILTER_LEVEL;


end

