function [DSER,IR_Late] = iceberg_core(IR)

    % Native equivalent of ita_split(IR,1) for analytical purposes
    omnichannelIR = IR;
    omnichannelIR.time = IR.time(:, 1);
    omnichannelIR.nChannels = 1;
    
    [IR_Early, shiftIndex] = native_time_shift(IR, 'auto'); % Pull FULL array
    
    try
        % Compute Center Time (Ts) using ISO 3382 native routine on mono reference
        cTime = native_center_time(omnichannelIR);
        
        IR_Early = native_time_window(IR_Early, [0 cTime], 'time', 'hann');
        % Shift back to fit the future composition
        DSER = native_time_shift(IR_Early, abs(shiftIndex));            % Push
    catch
        disp('ATTENTION: Can not calculate Center Time')
        IR_Early = native_time_window(IR_Early, [0 .01], 'time', 'rectwin'); 
        DSER = native_time_shift(IR_Early, abs(shiftIndex));            % Push
        cTime = 0.01;
    end   

    %% Shift IR to the begining
    [IR_Late, shiftIndex] = native_time_shift(IR, 'auto');
    
    % Crop out the Direc Sound+Early Reflections
    IR_Late = native_time_crop(IR_Late, [cTime 0], 'time');
    
    % Shift back to fit the future composition
    resyncSamples = IR.nSamples - IR_Late.nSamples;               
    IR_Late.time = [zeros(resyncSamples, IR_Late.nChannels); IR_Late.time];      % Adjusting time
    IR_Late.nSamples = size(IR_Late.time, 1);
    IR_Late.trackLength = IR_Late.nSamples / IR_Late.samplingRate;
    
    IR_Late = native_time_window(IR_Late, ...                       
        [0 (IR_Late.trackLength-0.05)], 'time', ...  
        'rectwin');                   
    IR_Late = native_time_shift(IR_Late, abs(shiftIndex));          % Push
end