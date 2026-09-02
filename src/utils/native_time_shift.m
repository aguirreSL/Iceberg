function [actVal, shiftAmountRet] = native_time_shift(signalObj, shiftAmount, varargin)
    % Native replacement for ita_time_shift
    %
    % Modes:
    %   native_time_shift(s, 'auto')              onset-shift (ita_start_IR ISO 3382, 20 dB)
    %   native_time_shift(s, '20dB')              onset-shift with explicit threshold
    %   native_time_shift(s, n, 'samples')        explicit circular shift by n samples
    %   native_time_shift(s, t, 'time')           explicit circular shift by t seconds
    %
    % The old heuristic that treated a numeric shift > 100 as samples and
    % anything else as seconds multiplied small sample shifts by fs and
    % corrupted them; callers must now pass the mode explicitly. 'auto' uses
    % the ISO 3382 onset (ita_time_shift/ita_start_IR semantics: earliest
    % channel onset across the signal), not the absolute peak.

    mode = 'time';
    if ~isempty(varargin) && ischar(varargin{1})
        mode = varargin{1};
    elseif ischar(shiftAmount)
        if strcmpi(shiftAmount, 'auto')
            mode = 'auto';
        elseif numel(shiftAmount) >= 2 && strcmpi(shiftAmount(end-1:end), 'dB')
            mode = 'auto';
        else
            error('NATIVE_TIME_SHIFT: unknown mode string %s', shiftAmount);
        end
    end

    fs = signalObj.samplingRate;
    timeData = signalObj.time;
    numCh = size(timeData, 2);

    % Auto / threshold shift: ISO 3382 onset across channels (ITA uses
    % min over channels, then shift_samples = -start + 1)
    if strcmpi(mode, 'auto')
        if numel(shiftAmount) >= 2 && strcmpi(shiftAmount(end-1:end), 'dB')
            threshold = str2double(shiftAmount(1:end-2));
        else
            threshold = 20;
        end
        foundStart = native_start_IR(signalObj, threshold);
        shiftSamples = -min(foundStart) + 1;
        shiftAmountRet = shiftSamples;
    else
        if strcmpi(mode, 'time')
            shiftSamples = round(shiftAmount * fs);
        elseif strcmpi(mode, 'samples')
            shiftSamples = shiftAmount;
            if shiftSamples ~= round(shiftSamples)
                warning('NATIVE_TIME_SHIFT: non-integer sample shift; rounding by circshift.');
            end
        else
            error('NATIVE_TIME_SHIFT: unknown mode %s', mode);
        end
        shiftAmountRet = shiftSamples;
    end

    % Perform circular shift exactly matching ITA algorithm
    shiftedData = zeros(size(timeData));
    for ch = 1:numCh
        shiftedData(:, ch) = circshift(timeData(:, ch), [double(shiftSamples) 0]);
    end

    actVal = signalObj;
    actVal.time = shiftedData;
    actVal.shiftIndex = shiftAmountRet;
end
