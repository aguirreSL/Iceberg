function calibrated_vbap = calibrate_vbap(signal_to_play, level, iAngles, configurationSetup)
% CALIBRATE_VBAP  Apply per-LS frequency filter and SPL-aligned level to a
% VBAP-rendered multi-channel signal. Native struct version.

iFs          = signal_to_play.samplingRate;
nLen         = signal_to_play.nSamples;
Level_Factor = configurationSetup.newLevelFactor;
lsdBperVolt  = 20 * log10(configurationSetup.iFactor / 2e-5);

%% Map ls_dir to physical LS indices in the master array
nLS = length(configurationSetup.ls_dir);
activeLSNumbers = zeros(1, nLS);
for ii = 1:nLS
    activeLSNumbers(ii) = find(configurationSetup.lsArray == configurationSetup.ls_dir(ii,1), 1);
end

%% Local frequency vector (replaces itaAudio's signal.freqVector)
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

%% Two nearest LS to the source angle.
% Replaces the hardcoded 0/90/180/270 cascade. The previous cascade had
% s1/s2 inverted in 3 of the 8 octants (sent higher SPL to the *further*
% LS); the defensive deal(max,min) on the levels did not fix it because it
% swapped values but not channel pointers.
allAngles = configurationSetup.ls_dir(:,1);
angDist   = abs(mod(allAngles - iAngles + 180, 360) - 180);
[~, order] = sort(angDist);
s1 = activeLSNumbers(order(1));
s2 = activeLSNumbers(order(2));

gap = abs(mod(allAngles(order(2)) - allAngles(order(1)) + 180, 360) - 180);
if gap == 0
    ratio = 0;
else
    ratio = angDist(order(1)) / gap;
end
s1_level = cos(ratio * pi/2)^2;
s2_level = sin(ratio * pi/2)^2;

%% Per-channel SPL targets
levels = zeros(1, max(activeLSNumbers));
if level ~= 'n'
    s1Db = 20 * log10(10^(level/20) * s1_level);
    s2Db = 20 * log10(10^(level/20) * s2_level);
    if ~isinf(s1Db), levels(s1) = s1Db; end
    if ~isinf(s2Db), levels(s2) = s2Db; end
end

%% RMS normalize, scale to target dB, apply per-LS spectral filter
out = zeros(nLen, length(activeLSNumbers));
for idx = 1:length(activeLSNumbers)
    ch = signal_to_play.time(:, idx);
    if any(isnan(ch))
        fprintf('NaN Channel, take a look ch: %i\n', idx);
        ch(isnan(ch)) = 0;
    end

    scaler = sqrt(mean(ch.^2)) * 2;
    if scaler ~= 0
        ch = ch / scaler;
    end

    ch = ch * 10^((levels(activeLSNumbers(idx)) - lsdBperVolt) / 20);

    chFFT = fft(ch);
    filterResp = Interpolation(:, activeLSNumbers(idx));
    if length(filterResp) ~= length(chFFT)
        filterResp = interp1( ...
            linspace(0, 1, length(filterResp)), filterResp, ...
            linspace(0, 1, length(chFFT))).';
    end
    out(:, idx) = real(ifft(chFFT .* filterResp)) * Level_Factor(activeLSNumbers(idx));
end

calibrated_vbap.time         = out;
calibrated_vbap.samplingRate = iFs;
calibrated_vbap.nChannels    = size(out, 2);
calibrated_vbap.nSamples     = nLen;
calibrated_vbap.trackLength  = nLen / iFs;
end
