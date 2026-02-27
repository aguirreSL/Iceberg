function actVal = native_time_window(signalObj, timeRange, mode, varargin)
    % Native replacement for ita_time_window
    
    if ~strcmp(mode, 'time')
        error('Only time mode supported currently in native replacement');
    end
    
    % Parse window type (default to hanning if rectangular not specified)
    windowType = 'hann';
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'windowType')
            windowType = varargin{i+1};
        end
    end
    
    fs = signalObj.samplingRate;
    
    startTime = timeRange(1);
    endTime = timeRange(2);
    
    startSampleRaw = round(startTime * fs) + 1;
    if endTime == 0
        endSampleRaw = signalObj.nSamples + 1;
    else
        endSampleRaw = round(endTime * fs) + 1;
    end
    
    % ITA's hidden logic: forces the start and end indices to be purely even numbers
    startSample = startSampleRaw + mod(startSampleRaw, 2);
    endSample = endSampleRaw + mod(endSampleRaw, 2);
    
    % PDI bug fix embedded in ITA source for signal lengths
    if endSample == signalObj.nSamples + 1
        endSample = endSample - 1;
    end
    
    startSample = max(1, min(startSample, signalObj.nSamples));
    endSample = max(1, min(endSample, signalObj.nSamples));
    
    windowLength = endSample - startSample + 1;
    
    if windowLength <= 0
        actVal = signalObj;
        actVal.time = [];
        actVal.nSamples = 0;
        actVal.trackLength = 0;
        return;
    end
    
    % ITA's standard time_window [start end] with 'time' mode keeps everything 
    % before 'start' normal, and fades out between 'start' and 'end' using 
    % the right half of a hanning window. Everything after 'end' is 0.
    
    if strcmpi(windowType, 'rectwin')
        win = ones(windowLength, 1);
    else
        % ITA uses window(@hann, lengthWindow) which calls hann(), NOT hanning().
        % CRITICAL: hann(N) and hanning(N) use DIFFERENT trig formulas:
        %   hann(N):    w(k) = 0.5*(1 - cos(2*pi*k/(N-1)))   (symmetric, endpoints = 0)
        %   hanning(N): w(k) = 0.5*(1 - cos(2*pi*k/(N+1)))   (periodic-like, no zero endpoints)
        % Using hanning() instead of hann() causes a 1.89e-4 systematic deviation.
        itaDefinedLength = 2 * (abs(endSample - startSample) + 1);
        fullOuterWin = hann(itaDefinedLength); 
        
        % Left side 1, right side fades to 0
        % Extracted exact array subset catching the tail "falling window slope"
        win = fullOuterWin(end - windowLength + 1 : end);
    end
    
    timeData = signalObj.time;
    
    % Apply fading window
    for ch = 1:size(timeData, 2)
        timeData(startSample:endSample, ch) = timeData(startSample:endSample, ch) .* win;
        
        % Keep everything before startSample UNTOUCHED (1.0 factor implicitly)
        
        % Setting purely outside of window end to 0
        timeData(endSample+1:end, ch) = 0;
    end
    
    actVal = signalObj;
    actVal.time = timeData;
end
