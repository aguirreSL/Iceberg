function actVal = native_time_crop(signalObj, timeRange, mode)
    % Native replacement for ita_time_crop, replicating its exact semantics:
    %  - interval = round([t1 t2]*fs) + 1
    %  - odd crop lengths are forced even (interval(2) decremented)
    %  - inverted ranges such as [t 0] do NOT keep t1..end; ITA keeps
    %    [1:(i2-1), (i1+1):end], i.e. removes the block in between (for
    %    [t 0] the result starts one sample later than a naive t1:end crop).

    if ~strcmp(mode, 'time')
        error('Only time mode supported currently in native replacement');
    end

    fs = signalObj.samplingRate;

    interval = round(timeRange .* fs) + 1;

    d = interval(2) - interval(1);
    if d > 0 && mod(d, 2) == 0
        interval(2) = interval(2) - 1;   % ITA: prevent odd sample numbers
    end

    if interval(2) > interval(1)
        timeData = signalObj.time(interval(1):interval(2), :);
    else
        timeData = signalObj.time([1:(interval(2)-1), (interval(1)+1):end], :);
    end

    actVal.time = timeData;
    actVal.samplingRate = fs;
    actVal.nSamples = size(timeData, 1);
    actVal.nChannels = size(timeData, 2);
    actVal.trackLength = actVal.nSamples / fs;
end
