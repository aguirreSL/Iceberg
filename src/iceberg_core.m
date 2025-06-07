function [DSER,IR_Late] = iceberg_core(IR)

omnichannelIR = ita_split(IR,1);                                %Select/get Omnidirectional IR
[IR_Early, shiftIndex] = ita_time_shift(omnichannelIR,'auto'); % Pull

try
    t20_bb = ita_roomacoustics(ita_normalize_dat(IR.ch(1)),'T20','broadbandanalysis',true);
    t20 = ita_roomacoustics(ita_normalize_dat(IR.ch(1)),'T20');
    sprintf('T20 Broadband %.3f s\nT20 Max %.3f s',t20_bb.T20.freqData,max(t20.T20.freqData))
    
    centerTime = ita_roomacoustics(omnichannelIR,'Center_Time','broadbandAnalysis',1);
    cTime = centerTime.Center_Time.freq;
    if isnan(cTime)
        centerTime = ita_roomacoustics(ita_normalize_dat(omnichannelIR),'Center_Time','broadbandAnalysis',1,'edcMethod','noCut');
        cTime = centerTime.Center_Time.freq;
    end
    IR_Early = ita_time_window(IR_Early,[0 cTime],'time','windowType','hann');
    % Shift back to fit the future composition
    DSER = ita_time_shift(IR_Early,abs(shiftIndex));            % Push
catch
    disp('ATTENTION: Can not calculate T20')
    IR_Early = ita_time_window(IR_Early,[0 .01],'time','windowType','rectwin'); %It does not make sense the cTime Only DS from Odeon simulation
    DSER = ita_time_shift(IR_Early,abs(shiftIndex));            % Push
    cTime = 0.01;
end   
    %% Shift IR to the begining
    [IR_Late, shiftIndex] = ita_time_shift(IR,'auto');
    % Crop out the the Direc Sound+Early Reflections (Center Time) (May work
    % with window as well, but the energy balance need to be verified in this case
    IR_Late = ita_time_crop(IR_Late,[cTime 0],'time');
    % Shift back to fit the future composition
    resyncSamples = IR.nSamples-IR_Late.nSamples;               % Timefactor
    IR_Late.time = [zeros(resyncSamples,4); IR_Late.time];      % Adjusting time
    IR_Late = ita_time_window(IR_Late,...                       %Get rid of non linearities
        [0 (IR_Late.trackLength-0.05)],'time',...  
        'windowType','rectwin');                   
    IR_Late = ita_time_shift(IR_Late,abs(shiftIndex));          % Push
end