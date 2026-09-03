function actVal = native_convolve(signalObj, irObj)
    % Native replacement for ita_convolve, replicating its algorithm choice:
    %  - fftDegree difference < 2  -> overlap-add via fftfilt (as ITA)
    %  - otherwise                 -> linear convolution by spectral
    %                                multiplication (as ITA's source*filter)
    % Output length is ITA's finalSamples = 2*ceil((n1+n2-1)/2), i.e. the
    % linear convolution length forced even. Channel pairing follows ITA's
    % fftfilt behaviour: channel k of the output pairs the k-th (or the only)
    % channel of each input.

    fs = signalObj.samplingRate;
    if fs ~= irObj.samplingRate
        error('Sample rates must match for convolution.');
    end

    numSignalCh = size(signalObj.time, 2);
    numIrCh = size(irObj.time, 2);

    finalSamples = size(signalObj.time, 1) + size(irObj.time, 1) - 1;
    finalSamples = 2 * ceil(finalSamples / 2);   % ITA: should be even for the fft

    overlapAdd = abs(log2(size(signalObj.time,1)) - log2(size(irObj.time,1))) < 2;

    timeData = zeros(finalSamples, max(numSignalCh, numIrCh));

    for k = 1:max(numSignalCh, numIrCh)
        x = zeros(finalSamples, 1);
        x(1:size(signalObj.time,1)) = signalObj.time(:, min(k, numSignalCh));
        h = zeros(finalSamples, 1);
        h(1:size(irObj.time,1)) = irObj.time(:, min(k, numIrCh));

        if overlapAdd
            convResult = fftfilt(h, x);
        else
            convResult = real(ifft(fft(x) .* fft(h)));
        end
        timeData(:, k) = convResult;
    end

    actVal.time = timeData;
    actVal.samplingRate = fs;
    actVal.nSamples = size(timeData, 1);
    actVal.nChannels = size(timeData, 2);
    actVal.trackLength = actVal.nSamples / fs;
end
