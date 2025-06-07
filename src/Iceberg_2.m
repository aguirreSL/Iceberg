function [final_audio] = Iceberg(signal,level,angle,iIR,configurationSetup)
%% Iceberg2021(signal,level,angle,iIR,configurationSetup)
% Auralize file to a specific array in hybrid mode VBAP/Ambisonics
% Three IR Options are provided
 

% IR's Path
if ispc
    IR00Path = [pwd '\wavFiles\rt00\Restaurant.BFormat'];
    IR05Path = [pwd '\wavFiles\rt05\rum019.BFormat'];
    IR11Path = [pwd '\wavFiles\rt11\Restaurant.BFormat'];
elseif ismac
    IR00Path = [pwd '/wavFiles/rt00/Restaurant.BFormat'];
    IR05Path = [pwd '/wavFiles/rt05/rum019.BFormat'];
    IR11Path = [pwd '/wavFiles/rt11/Restaurant.BFormat'];

end
%% initialize itaAudio objects
VBAP_DSER_Part = itaAudio();
Amb_ERLR_Part = itaAudio();
VBAP_DSER_Part.samplingRate = signal.samplingRate;
Amb_ERLR_Part.samplingRate = signal.samplingRate;
%% Array Settings
activeLSNumbers = zeros(1,length(configurationSetup.ls_dir));
for indexx= 1:length(configurationSetup.ls_dir)
    iArrayP = find(configurationSetup.LSArray==configurationSetup.ls_dir(indexx,1),1);
    activeLSNumbers(indexx) = iArrayP;
end

iFs         = signal.samplingRate;                                      %Sample Frequency
signal      = ita_normalize_dat(signal,'allchannels','true');
signalAmb   = set_level_ambisonics_fly_in(signal,level,angle,configurationSetup);    %Normalizing all
[D4,~] = ambiDecoder(configurationSetup.ls_dir,'SAD',1,1);                                 % Create Decoder Matrix with n=4 LS SAD decoder! 
                                                                        %   - Sampling decoder (SAD)
                                                                        %   - Mode-matching decoder (MMD)
                                                                        %   - Energy-preserving decoder (EPAD)
                                                                        %   - All-round ambisonic panning (ALLRAD)
                                                                        %   - Constant Angular Spread Decoder (CSAD)
%% Load the simulated RIR
iOdeonPosition = (angle/5)+1;                                  %Specific angles from Odeon Simulation
if iIR == 1
    IR = ita_read([IR05Path num2str(iOdeonPosition) '.Wav']);   %Load Ambisonics IR 0.5s
elseif iIR ==2
    IR = ita_read([IR11Path num2str(iOdeonPosition) '.Wav']);   %Load Ambisonics IR 1.1s
elseif iIR ==3
    IR = ita_read([IR00Path num2str(iOdeonPosition) '.Wav']);   %Load Ambisonics IR anechoic
end

omnichannelIR = ita_split(IR,1);                                %Select/get Omnidirectional IR

%% Center time 
[IR_Early, shiftIndex] = ita_time_shift(omnichannelIR,'auto'); % Pull

if iIR == 3
    IR_Early = ita_time_window(IR_Early,[0 .01],'time','windowType','rectwin'); %It does not make sense the cTime Only DS from Odeon simulation
    DSER = ita_time_shift(IR_Early,abs(shiftIndex));            % Push
    cTime = 0.01;
else
    centerTime = ita_roomacoustics(omnichannelIR,'Center_Time','broadbandAnalysis',1);
    cTime = centerTime.Center_Time.freq;
    if isnan(cTime)
        centerTime = ita_roomacoustics(ita_normalize_dat(omnichannelIR),'Center_Time','broadbandAnalysis',1,'edcMethod','noCut');
        cTime = centerTime.Center_Time.freq;
    end
    IR_Early = ita_time_window(IR_Early,[0 cTime],'time','windowType','hann');
    % Shift back to fit the future composition
    DSER = ita_time_shift(IR_Early,abs(shiftIndex));            % Push
end
    
    %% Signal with Omnidirectional reverberation
    convolved_DSER_signal = ita_convolve(signal,DSER);
    %% Calibrated VBAP to specified Sound Pressure Level
    VBAP_DS =   set_level_vbap_fly_in(convolved_DSER_signal,level,angle,configurationSetup);
    VBAP_DS =   VBAP_DS*max(DSER.time);
    %% fetch VBAP Direct Sound to the Array
    
    for i = 1:length(activeLSNumbers)
        VBAP_DSER_Part.time(:,activeLSNumbers(i)) = VBAP_DS.time(:,i);
    end
    %% Shift IR to the begining
    [IR_Late, shiftIndex] = ita_time_shift(IR,'auto');
    % Crop out the the Direc Sound+Early Reflections (Center Time) (May work
    % with window as well, but the energy balance need to be verified in this
    % case
    IR_Late = ita_time_crop(IR_Late,[cTime 0],'time');
    % Shift back to fit the future composition
    resyncSamples = IR.nSamples-IR_Late.nSamples;               % Timefactor
    IR_Late.time = [zeros(resyncSamples,4); IR_Late.time];      % Adjusting time
    IR_Late = ita_time_window(IR_Late,...                       %Get rid of non linearities
        [0 (IR_Late.trackLength-0.05)],'time',...  %Get rid of non linearities
        'windowType','rectwin');                   %Get rid of non linearities
    IR_Late = ita_time_shift(IR_Late,abs(shiftIndex));          % Push
    
    %% Ambisonics
    sinal_Ambisonics_ERLR   = ita_convolve(signalAmb,IR_Late);         % Signal with Late Reverberation
    %% Calibrated Ambisonics to specified Sound Pressure Level
    %Decode to Loudspeaker array
    Ambisonics_ERLRSignal = itaAudio(decodeBformat(...
        sinal_Ambisonics_ERLR.time,...
        D4),iFs,'time');
    % fetch Ambisonics Late reverberation to the Array
    for i = 1:length(activeLSNumbers)
        Amb_ERLR_Part.time(:,activeLSNumbers(i))  = Ambisonics_ERLRSignal.time(:,i);
    end
    %% Combine DS ER and LR
    final_audio = ita_add(VBAP_DSER_Part,Amb_ERLR_Part);    %%
    if length(configurationSetup.LSArray) >final_audio.dimensions
        final_audio.time(:,final_audio.dimensions+1:length(configurationSetup.LSArray)) = zeros(final_audio.nSamples,length(configurationSetup.LSArray)-(final_audio.dimensions));
    end
end
