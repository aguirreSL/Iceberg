function [actVal, shiftAmountRet] = native_time_shift(signalObj, shiftAmount, varargin)
    % Native replacement for ita_time_shift
    
    mode = 'time';
    if ~isempty(varargin) && ischar(varargin{1})
        mode = varargin{1};
    elseif isnumeric(shiftAmount) && shiftAmount > 100
        % Fallback intuition if missing explicit 'time' char
        mode = 'samples'; 
    end
    
    fs = signalObj.samplingRate;
    timeData = signalObj.time;
    numSamples = size(timeData, 1);
    numCh = size(timeData, 2);
    
    % Auto shift logic
    if ischar(shiftAmount) && strcmpi(shiftAmount, 'auto')
        % Simple auto shift: find max peak and shift it to the beginning or specifically as needed
        % Actually, auto shift in ITA usually aligns the impulse response peak to 0. 
        % We will return the shift amount to allow pushing back later
        [~, maxIdx] = max(abs(timeData(:, 1)));
        shiftSamples = -(maxIdx - 1); % Shift peak to front
        shiftAmountRet = shiftSamples; % Usually returned as shiftIndex
    else
        % It's a specific amount
        if strcmpi(mode, 'time')
            shiftSamples = round(shiftAmount * fs);
        else
            shiftSamples = shiftAmount;
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
    actVal.shiftIndex = shiftAmountRet; % Custom field to mimic multiple returns if needed
end
