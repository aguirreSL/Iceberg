function cTime = native_center_time(IR)
% NATIVE_CENTER_TIME Calculates the center time (Ts) of an impulse response
% mirroring the exact analytical structure of ita_roomacoustics_EDC.m.
%
% IR: audioStruct with .time and .samplingRate (must be single channel)

    fs = IR.samplingRate;
    p = IR.time;
    
    if size(p, 2) > 1
        p = p(:, 1); % force single channel for computation
    end
    
    if isrow(p)
        p = p';
    end
    
    nSamples = length(p);
    energyData = p.^2;
    timeVector = (0:nSamples-1)' / fs;
    
    % Smooth data for intersection finding (ITA uses 75ms blocks)
    smoothBlockLength = 0.075;
    nSamplesPerBlock = round(smoothBlockLength * fs);
    numBlocks = floor(nSamples / nSamplesPerBlock);
    
    if numBlocks < 2
        % Fallback for extremely short IR
        cTime = sum(energyData .* timeVector) / (sum(energyData) + eps);
        return;
    end
    
    energyTrunc = energyData(1:numBlocks*nSamplesPerBlock);
    timeWinData = sum(reshape(energyTrunc, nSamplesPerBlock, numBlocks), 1)' / nSamplesPerBlock;
    timeVecWin  = (0.5 + (0:numBlocks-1))' * nSamplesPerBlock / fs;
    
    % Estimate noise from the last 10%
    noise_idx = max(1, floor(0.9 * numBlocks));
    pSquareAtIntersection = mean(timeWinData(noise_idx:end));
    
    % Find intersection t1 idx
    t1idx = find(timeWinData > pSquareAtIntersection * 2, 1, 'last');
    if isempty(t1idx)
        t1idx = numBlocks;
    end
    
    t1 = timeVecWin(t1idx);
    t1IdxRaw = min(nSamples, round(t1 * fs) + 1);
    
    % Find t0 idx (10 dB above noise)
    t0idx = find(timeWinData(1:t1idx) > 10 * pSquareAtIntersection, 1, 'last');
    if isempty(t0idx)
        t0idx = max(1, t1idx - 5);
    end
    
    % Regression on smoothed data for TofLast10dB
    if t0idx < t1idx
        X = [timeVecWin(t0idx:t1idx).^0 timeVecWin(t0idx:t1idx)];
        coeff = X \ (10*log10(abs(timeWinData(t0idx:t1idx)) + eps));
        TofLast10dB = -60 / coeff(2);
    else
        TofLast10dB = 0;
    end
    
    if TofLast10dB < 0 || TofLast10dB > 20 || isnan(TofLast10dB)
        TofLast10dB = 0;
    end
    
    % Calculate correction C according to DIN EN ISO 3382
    C = pSquareAtIntersection * TofLast10dB / (6 * log(10)) * fs; 
    
    % Subtract Noise Energy
    energyDataClean = energyData - pSquareAtIntersection;
    
    % Negative bounding for EDC
    EDC = cumsum(energyDataClean(t1IdxRaw:-1:1));
    EDC = EDC(end:-1:1) + C;
    
    % Exact center time formula from ita_roomacoustics_EDC line 307
    numerator  = sum(energyDataClean(1:t1IdxRaw) .* timeVector(1:t1IdxRaw)) + C^2 + C * pSquareAtIntersection * t1;
    cTime = numerator / (EDC(1) + eps);
end
