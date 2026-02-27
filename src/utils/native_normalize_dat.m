function actVal = native_normalize_dat(signalObj)
    % Native replacement for ita_normalize_dat
    actVal = signalObj;
    
    if isempty(signalObj.time)
        return;
    end
    
    maxVal = max(abs(signalObj.time(:)));
    
    if maxVal > 0
        actVal.time = signalObj.time / maxVal;
    end
end
