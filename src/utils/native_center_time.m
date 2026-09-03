function cTime = native_center_time(IR)
% NATIVE_CENTER_TIME Centre time (Ts) of an impulse response, replicating the
% ITA-Toolbox default path bit for bit:
%   ita_roomacoustics preprocessing (ita_time_shift '20dB', wrapped-tail
%   zeroing, trailing-zero truncation), then the broadband Lundeby
%   estimation (ita_roomacoustics_reverberation_time_lundeby with 30 ms
%   windows), then the EDC centre-time formula of ita_roomacoustics_EDC with
%   method 'cutWithCorrection' (default).
%
% On Lundeby failure (NaN, as in ITA) the function falls back to the full
% uncut centroid, mirroring the thesis-era isnan(cTime) retry with
% 'edcMethod','noCut' (normalisation is irrelevant: the centroid is
% scale-invariant).
%
% IR: audioStruct with .time and .samplingRate (must be single channel)

    fs = IR.samplingRate;
    p = IR.time;
    if size(p, 2) > 1
        p = p(:, 1);
    end
    if isrow(p)
        p = p';
    end
    nSamples = numel(p);

    %% ---- ita_roomacoustics preprocessing ----
    % onset shift (ita_time_shift '20dB')
    shiftSamples = -(native_start_IR(IR, 20)) + 1;      % ITA: -start + 1
    p = circshift(p, [double(shiftSamples), 0]);

    % nInputSamples: original onset position in the shifted frame, even-forced
    nInputSamples = nSamples + shiftSamples + rem(shiftSamples, 2);

    % zero the wrapped-around tail, then walk back over trailing zeros
    tmpData = p;
    tmpData(nInputSamples:end) = 0;
    if nInputSamples > 1
        startValue = nInputSamples;
        while tmpData(nInputSamples - 1) == 0
            nInputSamples = nInputSamples - 1;
            if nInputSamples == 1
                error('native_center_time: channel is empty');
            end
        end
        if nInputSamples < 5
            nInputSamples = startValue;   % ITA: cancel cutting if < 5 samples left
        end
    end
    p = tmpData(1:nInputSamples);
    nSamples = numel(p);

    %% ---- broadband Lundeby (verbatim structure) ----
    [revT, noiseEst, intersectionTime] = native_lundeby_broadband(p, fs);

    %% ---- ita_roomacoustics_EDC 'cutWithCorrection' + calcCenterTime ----
    energyData = p.^2;
    timeVector = (0:nSamples-1)' / fs;

    if isnan(revT) || isnan(intersectionTime)
        % Lundeby failed; ITA returns NaN, the thesis chain retried with
        % 'edcMethod','noCut' = full uncut centroid (C = 0, t1IdxRaw = end).
        cTime = sum(energyData .* timeVector) / sum(energyData);
        return;
    end

    t1 = intersectionTime;
    [~, t1IdxRaw] = min(abs(timeVector - t1));

    pSquareAtIntersection = noiseEst;         % noiseRMS.^2 (linear energy)
    TofLast10dB = revT;
    C = pSquareAtIntersection * TofLast10dB / (6 * log(10)) * fs;

    EDCfirst = sum(energyData(1:t1IdxRaw)) + C;
    numerator = sum(energyData(1:t1IdxRaw) .* timeVector(1:t1IdxRaw)) ...
        + C^2 + C * pSquareAtIntersection * t1;
    cTime = numerator / EDCfirst;
end

function [revT, noiseEst, intersectionTime] = native_lundeby_broadband(rawData, fs)
    % Broadband branch of ita_roomacoustics_reverberation_time_lundeby:
    % 30 ms windows, iterative noise/decay/crossing estimation. Includes the
    % function's own onset shift (ita_time_shift '20dB' at its top, a no-op
    % for already-shifted data in practice) and its nSamples2cut trimming;
    % the returned intersectionTime is referenced to the input frame
    % (crossingPoint - timeShifted), as in ITA's output object.

    noiseEst = NaN;
    revT = NaN;
    intersectionTime = NaN;

    % internal onset shift (usually a no-op: input is already shifted)
    innerStart = native_start_IR(struct('time', rawData, 'samplingRate', fs, ...
        'nSamples', numel(rawData), 'nChannels', 1), 20);
    innerShift = -(innerStart(1)) + 1;
    timeShifted = innerShift / fs;
    rawData = circshift(rawData, [double(innerShift), 0]);

    nSamples2cut = -timeShifted * fs + rem(-timeShifted * fs, 2);
    nSamples = numel(rawData) - nSamples2cut;

    freqDepWinTime = 0.03;                            % broadband window size
    nPartsPer10dB = 5;
    dbAboveNoise = 10;
    useDynRangeForRegression = 20;

    rawTimeData = rawData.^2;

    % 1) smooth
    nSamplesPerBlock = round(freqDepWinTime * fs);
    timeWinData = local_blockMean(rawTimeData(1:nSamples), nSamplesPerBlock);
    timeVecWin = (0:size(timeWinData,1)-1)' * nSamplesPerBlock / fs;

    % 2) estimate noise
    noiseEst = mean(timeWinData(end-round(size(timeWinData,1)/10):end)) + realmin;

    % 3) regression
    [~, startIdx] = max(timeWinData);
    stopIdx = find(10*log10(timeWinData(startIdx+1:end)) > 10*log10(noiseEst) + dbAboveNoise, 1, 'last') + startIdx;
    if isempty(stopIdx)
        return
    end
    dynRange = diff(10*log10(timeWinData([startIdx stopIdx])));
    if (stopIdx == startIdx) || dynRange > -5
        return
    end

    X = [ones(stopIdx-startIdx+1,1) timeVecWin(startIdx:stopIdx)];
    c = X \ (10*log10(timeWinData(startIdx:stopIdx)));
    if c(2) == 0 || any(isnan(c))
        return
    end

    % 4) preliminary crossing point
    crossingPoint = (10*log10(noiseEst) - c(1)) / c(2);
    if crossingPoint > (numel(rawData) / fs + timeShifted) * 2
        return
    end

    % 5) new local time interval length
    nBlocksInDecay = diff(10*log10(timeWinData([startIdx stopIdx]))) / -10 * nPartsPer10dB;
    nSamplesPerBlock = round(diff(timeVecWin([startIdx stopIdx])) / nBlocksInDecay * fs);

    % 6) average
    timeWinData = local_blockMean(rawTimeData(1:nSamples), nSamplesPerBlock);
    timeVecWin = (0:size(timeWinData,1)-1)' * nSamplesPerBlock / fs;
    [~, idxMax] = max(timeWinData);

    oldCrossingPoint = 11 + crossingPoint;
    loopCounter = 0;
    while abs(oldCrossingPoint - crossingPoint) > 0.01
        % 7) estimate background level
        correspondingDecay = 10;
        idxLast10percent = round(size(timeWinData,1) * 0.9);
        idx10dBBelowCrosspoint = max(1, round((crossingPoint - correspondingDecay ./ c(2)) * fs / nSamplesPerBlock));
        noiseEst = mean(timeWinData(min(idxLast10percent, idx10dBBelowCrosspoint):end)) + realmin;

        % 8) estimate late decay slope
        startIdx = find(10*log10(timeWinData(idxMax:end)) < 10*log10(noiseEst) + dbAboveNoise + useDynRangeForRegression, 1, 'first') + idxMax - 1;
        if isempty(startIdx)
            startIdx = 1;
        end
        stopIdx = find(10*log10(timeWinData(startIdx+1:end)) < 10*log10(noiseEst) + dbAboveNoise, 1, 'first') + startIdx;
        if isempty(stopIdx)
            return
        end
        X = [ones(stopIdx-startIdx+1,1) timeVecWin(startIdx:stopIdx)];
        c = X \ (10*log10(timeWinData(startIdx:stopIdx)));
        if c(2) >= 0
            c(2) = Inf;
            break
        end

        % 9) find crosspoint
        oldCrossingPoint = crossingPoint;
        crossingPoint = (10*log10(noiseEst) - c(1)) / c(2);

        loopCounter = loopCounter + 1;
        if loopCounter > 30
            break
        end
    end

    revT = -60 / c(2);
    intersectionTime = crossingPoint - timeShifted;
end

function bm = local_blockMean(x, nSamplesPerBlock)
    % sum(reshape(x(1:floor(end/b)*b)), b, nBlocks), 1).' / b
    nBlocks = floor(numel(x) / nSamplesPerBlock);
    bm = sum(reshape(x(1:nBlocks*nSamplesPerBlock), nSamplesPerBlock, nBlocks), 1).' / nSamplesPerBlock;
end
