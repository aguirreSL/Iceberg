function actVal = native_add(signalObj1, signalObj2)
    % Native replacement for ita_add
    
    fs1 = signalObj1.samplingRate;
    fs2 = signalObj2.samplingRate;
    
    if fs1 ~= fs2
        error('Sample rates must match for addition.');
    end
    
    timeData1 = signalObj1.time;
    timeData2 = signalObj2.time;
    
    len1 = size(timeData1, 1);
    len2 = size(timeData2, 1);
    
    ch1 = size(timeData1, 2);
    ch2 = size(timeData2, 2);
    
    maxLen = max(len1, len2);
    maxCh = max(ch1, ch2);
    
    % Pad matrices to match size
    paddedData1 = zeros(maxLen, maxCh);
    paddedData1(1:len1, 1:ch1) = timeData1;
    
    paddedData2 = zeros(maxLen, maxCh);
    paddedData2(1:len2, 1:ch2) = timeData2;
    
    addedData = paddedData1 + paddedData2;
    
    actVal = signalObj1;
    actVal.time = addedData;
    actVal.nSamples = maxLen;
    actVal.nChannels = maxCh;
    actVal.trackLength = maxLen / fs1;
end
