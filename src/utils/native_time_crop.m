function actVal = native_time_crop(signalObj, timeRange, mode)
    % Native replacement for ita_time_crop
    
    if ~strcmp(mode, 'time')
        error('Only time mode supported currently in native replacement');
    end
    
    fs = signalObj.samplingRate;
    
    % timeRange is [startEndInSeconds]
    startSample = round(timeRange(1) * fs) + 1;
    if timeRange(2) == 0
        endSample = signalObj.nSamples;
    else
        endSample = round(timeRange(2) * fs);
    end
    
    % Prevent index bounds error
    startSample = max(1, min(startSample, signalObj.nSamples));
    endSample = max(1, min(endSample, signalObj.nSamples));
    
    if startSample > endSample
        timeData = [];
    else
        timeData = signalObj.time(startSample:endSample, :);
    end
    
    actVal.time = timeData;
    actVal.samplingRate = fs;
    actVal.nSamples = size(timeData, 1);
    actVal.nChannels = size(timeData, 2);
    actVal.trackLength = actVal.nSamples / fs;
end
