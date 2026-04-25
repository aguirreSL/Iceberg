function [DSER,IR_Late] = iceberg_core(IR)

    % Mono omnidirectional IR (W channel) — used both for cTime and as the
    % source for DSER. This matches the original ITA pipeline (ita_split(IR,1));
    % VBAP needs an envelope-only mono signal, not the full B-Format.
    omnichannelIR = IR;
    omnichannelIR.time = IR.time(:, 1);
    omnichannelIR.nChannels = 1;

    [IR_Early, shiftIndex] = native_time_shift(omnichannelIR, 'auto');

    try
        cTime = native_center_time(omnichannelIR);
        IR_Early = native_time_window(IR_Early, [0 cTime], 'time', 'windowType', 'hann');
        DSER = native_time_shift(IR_Early, abs(shiftIndex));
    catch
        disp('ATTENTION: Can not calculate Center Time')
        IR_Early = native_time_window(IR_Early, [0 .01], 'time', 'windowType', 'rectwin');
        DSER = native_time_shift(IR_Early, abs(shiftIndex));
        cTime = 0.01;
    end

    %% Late part keeps the full B-Format (4 ch) for Ambisonics decoding
    [IR_Late, shiftIndex] = native_time_shift(IR, 'auto');
    IR_Late = native_time_crop(IR_Late, [cTime 0], 'time');

    resyncSamples = IR.nSamples - IR_Late.nSamples;
    IR_Late.time = [zeros(resyncSamples, IR_Late.nChannels); IR_Late.time];
    IR_Late.nSamples = size(IR_Late.time, 1);
    IR_Late.trackLength = IR_Late.nSamples / IR_Late.samplingRate;

    IR_Late = native_time_window(IR_Late, ...
        [0 (IR_Late.trackLength-0.05)], 'time', ...
        'windowType', 'rectwin');
    IR_Late = native_time_shift(IR_Late, abs(shiftIndex));
end