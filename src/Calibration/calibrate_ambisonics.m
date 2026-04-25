function signalcalibrated = calibrate_ambisonics(signal, level, iAngles, configurationSetup)
% CALIBRATE_AMBISONICS  Apply per-LS frequency filter and SPL alignment to an
% Ambisonics signal using the nearest physical LS as virtual reference.
% Native struct version.

if level > 95
    error('Level max is 95 dB');
end

iFs         = signal.samplingRate;
nLen        = signal.nSamples;
lsdBperVolt = 20 * log10(configurationSetup.iFactor / 2e-5);

%% Map ls_dir to physical LS indices in the master array
nLS = length(configurationSetup.ls_dir);
activeLSNumbers = zeros(1, nLS);
for ii = 1:nLS
    activeLSNumbers(ii) = find(configurationSetup.lsArray == configurationSetup.ls_dir(ii,1), 1);
end

%% Local frequency vector
nFFT  = nLen;
nBins = floor(nFFT/2) + 1;
freqVec = (0:nBins-1)' * (iFs / nFFT);

Interpolation = zeros(nBins, max(activeLSNumbers));
for iCount = activeLSNumbers
    Interpolation(:, iCount) = pchip( ...
        configurationSetup.iLoudspeakerFreqFilter(iCount).freqVector, ...
        configurationSetup.iLoudspeakerFreqFilter(iCount).freq, ...
        freqVec);
end

%% Nearest physical LS as virtual reference (Nearest Speaker Pan).
% Replaces the hardcoded 0/90/180/270 cascade — works for any layout.
allAngles = configurationSetup.ls_dir(:,1);
angDist   = abs(mod(allAngles - iAngles + 180, 360) - 180);
[~, order] = sort(angDist);
iChannel = activeLSNumbers(order(1));

%% RMS normalize, scale to target dB, apply nearest-LS spectral filter
ch = signal.time;
ch(isnan(ch)) = 0;

scaler = sqrt(mean(ch.^2)) * 2;
if scaler ~= 0
    ch = ch / scaler;
end

if level ~= 'n'
    ch = ch * 10^((level - lsdBperVolt) / 20);
end

chFFT = fft(ch);
filterResp = Interpolation(:, iChannel);
if length(filterResp) ~= length(chFFT)
    filterResp = interp1( ...
        linspace(0, 1, length(filterResp)), filterResp, ...
        linspace(0, 1, length(chFFT))).';
end
filtered = real(ifft(chFFT .* filterResp)) * configurationSetup.newLevelFactor(iChannel);

signalcalibrated.time         = filtered;
signalcalibrated.samplingRate = iFs;
signalcalibrated.nChannels    = 1;
signalcalibrated.nSamples     = nLen;
signalcalibrated.trackLength  = nLen / iFs;
end
