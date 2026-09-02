function sampleStart = native_start_IR(IR, threshold)
    % NATIVE_START_IR Find the start of an impulse response (ISO 3382 A.3.4)
    % Native replica of ita_start_IR (ITA-Toolbox), local_search_ISO3382 path.
    %
    % IR: audioStruct with .time and .samplingRate (any channel count)
    % threshold: level below the maximum (dB, energy domain), default 20
    %
    % Returns a 1 x nChannels vector of onset sample indices. ita_time_shift
    % 'auto'/'XdB' shifts by -(min(sampleStart)) + 1 across channels.

    if nargin < 2
        threshold = 20;
    end

    timeData = IR.time;
    nSamples = size(timeData, 1);
    nCh      = size(timeData, 2);

    IRsquare = timeData.^2;

    % assume the last 10% of the IR is noise, and calculate its noise level
    NoiseLevel = mean(IRsquare(round(.9*nSamples):end, :), 1);

    % get the maximum of the signal, that is the assumed IR peak
    [max_val, max_idx] = max(IRsquare, [], 1);

    % bad-SNR guard: less than 20 dB SNR or peak in the noisy part -> no shift
    idxNoShift = max_val < 100*NoiseLevel | max_idx > round(0.9*nSamples);

    threshold = -abs(threshold);
    sampleStart = ones(1, nCh);

    for idx = 1:nCh
        if idxNoShift(idx)
            continue
        end

        % if maximum lies on the first point, there is nothing to search
        if max_idx(idx) > 1
            abs_dat = 10*log10(IRsquare(1:max_idx(idx), idx)) - 10*log10(max_val(idx));

            lastBelowThreshold = find(abs_dat < threshold, 1, 'last');
            if ~isempty(lastBelowThreshold)
                sampleStart(idx) = lastBelowThreshold;
            else
                sampleStart(idx) = 1;
            end

            % Check if oscillations exist before the last value below
            % threshold. If so, these are part of the RIR and need to be
            % considered.
            idx6dBaboveThreshold = find(abs_dat(1:sampleStart(idx)) > threshold + 6);
            if ~isempty(idx6dBaboveThreshold)
                tmp = find(abs_dat(1:idx6dBaboveThreshold(1)) < threshold, 1, 'last');
                if isempty(tmp)
                    sampleStart(idx) = 1;
                else
                    sampleStart(idx) = tmp;
                end
            end
        end
    end
end
