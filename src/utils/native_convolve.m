function actVal = native_convolve(signalObj, irObj)
    % Native replacement for ita_convolve
    % Note: Assumes both inputs are standard audio objects matching ITA structures
    
    fs = signalObj.samplingRate;
    if fs ~= irObj.samplingRate
        error('Sample rates must match for convolution.');
    end
    
    numSignalCh = size(signalObj.time, 2);
    numIrCh = size(irObj.time, 2);
    
    % If doing full matrix convolution
    timeData = zeros(size(signalObj.time, 1) + size(irObj.time, 1), max(numSignalCh, numIrCh));
    
    for k = 1:max(numSignalCh, numIrCh)
        sIdx = min(k, numSignalCh);
        iIdx = min(k, numIrCh);
        
        % Using built-in conv
        convResult = conv(signalObj.time(:, sIdx), irObj.time(:, iIdx));
        timeData(1:length(convResult), k) = convResult;
    end
    
    % Structure similar to itaAudio
    actVal.time = timeData;
    actVal.samplingRate = fs;
    actVal.nSamples = size(timeData, 1);
    actVal.nChannels = size(timeData, 2);
    actVal.trackLength = actVal.nSamples / fs; 
end
