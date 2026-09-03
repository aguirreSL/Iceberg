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

%% Pair selection and in-pair levels.
% configSetup.pairSelection (optional, default 'cascade2022'):
%   'cascade2022' - bit-faithful replication of the thesis-era
%       set_level_vbap_fly_in octant cascade, INCLUDING its inversion in
%       three of the eight octants (for sources in (90,180] and (225,270]
%       the louder coefficient goes to the FARTHER loudspeaker of the pair).
%       This is the parity target of the migration: reproduce the measured
%       chain exactly, error included, until the planned redesign removes
%       this stage altogether. Requires the 4-LS cardinal layout.
%   'nearest'     - corrected nearest-two selection (25 Apr 2026 fix).
% The in-pair level law itself (cos^2/sin^2 of the fractional position
% inside the pair, max to s1) is identical in both modes.
allAngles = configurationSetup.ls_dir(:,1);
angDist   = abs(mod(allAngles - iAngles + 180, 360) - 180);
[~, order] = sort(angDist);

if isfield(configurationSetup, 'pairSelection') && strcmpi(configurationSetup.pairSelection, 'nearest')
    s1 = activeLSNumbers(order(1));
    s2 = activeLSNumbers(order(2));
else
    % Verbatim octant cascade of set_level_vbap_fly_in.m (2022). Note the
    % (180,270] octant assigns s1=180/s2=270 in BOTH halves - that is the
    % original code, not a typo.
    lsOf = @(a) activeLSNumbers(find(allAngles == a, 1));
    cardinals = [0 90 180 270];
    if any(arrayfun(@(a) isempty(find(allAngles == a, 1)), cardinals))
        error('calibrate_vbap:cascade2022 requires the 4-LS cardinal layout (0/90/180/270).');
    end
    if iAngles > 90 && iAngles <= 180
        if iAngles <= 135
            s1 = lsOf(180); s2 = lsOf(90);
        else
            s1 = lsOf(90);  s2 = lsOf(180);
        end
    elseif iAngles > 180 && iAngles <= 270
        s1 = lsOf(180); s2 = lsOf(270);
    elseif iAngles > 270 && iAngles <= 360
        if iAngles <= 315
            s1 = lsOf(270); s2 = lsOf(0);
        else
            s1 = lsOf(0);   s2 = lsOf(270);
        end
    else
        if iAngles <= 45
            s1 = lsOf(0);   s2 = lsOf(90);
        else
            s1 = lsOf(90);  s2 = lsOf(0);
        end
    end
end

gap = abs(mod(allAngles(order(2)) - allAngles(order(1)) + 180, 360) - 180);
if gap == 0
    ratio = 0;
else
    ratio = angDist(order(1)) / gap;
end
s1_level = cos(ratio * pi/2)^2;
s2_level = sin(ratio * pi/2)^2;
% The 2022 code applied deal(max, min): the louder coefficient always goes
% to the cascade's s1 (in 'nearest' mode s1 is the nearest, so this is a
% no-op there).
lv = max(s1_level, s2_level);
lq = min(s1_level, s2_level);

%% Per-channel SPL targets
levels = zeros(1, max(activeLSNumbers));
if level ~= 'n'
    s1Db = 20 * log10(10^(level/20) * lv);
    s2Db = 20 * log10(10^(level/20) * lq);
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
    ch = double(ch);
    % lsdBperVolt is single when iFactor is (the .mat stores it so), which
    % makes the scaled stimulus single. itaAudio casts time data to double
    % before its fft; replicate that so the FFT runs in double, as in ITA.

    chFFT = fft(ch);
    filterResp = double(Interpolation(:, activeLSNumbers(idx)));
    % Cast to double before the multiply, as itaAudio does when constructed
    % from (possibly single-precision) freq data: the calibration .mat
    % stores the EQ curves as single, and multiplying a double FFT by a
    % single filter collapses the whole product to single precision.
    % Mirror the half-spectrum onto the full FFT grid (Hermitian symmetry),
    % as ita_multiply_spk does. The previous code stretched the nBins filter
    % over the nFFT bins, replacing H(f) with 0.5*(H(f/2)+H(fs/2-f/2)) and
    % smearing steep curve sections by up to ~4.7 dB per band.
    nFFTlen = numel(chFFT);
    nBinsLen = numel(filterResp);
    if nBinsLen ~= nFFTlen
        if mod(nFFTlen, 2) == 0
            filterResp = [filterResp; conj(flipud(filterResp(2:end-1)))];
        else
            filterResp = [filterResp; conj(flipud(filterResp(2:end)))];
        end
    end
    out(:, idx) = real(ifft(chFFT .* filterResp)) * Level_Factor(activeLSNumbers(idx));
end

calibrated_vbap.time         = out;
calibrated_vbap.samplingRate = iFs;
calibrated_vbap.nChannels    = size(out, 2);
calibrated_vbap.nSamples     = nLen;
calibrated_vbap.trackLength  = nLen / iFs;
end
