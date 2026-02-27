function calibrated_vbap = calibrate_vbap(signal_to_play, configurationSetup)
%% Load values from the most recent calibration
Level_Factor = configurationSetup.newLevelFactor; %LS Level adjust
iFs = signal.samplingRate;                          %Sample Frequency
lsdBperVolt = (20*log10((configurationSetup.iFactor)/2e-5)); %Factor

activeLSNumbers = zeros(1,length(configurationSetup.ls_dir));
for indexx= 1:length(configurationSetup.ls_dir)
iArrayP = find(configurationSetup.LSArray==configurationSetup.ls_dir(indexx,1),1);
activeLSNumbers(indexx) = iArrayP;
end
for iCount = activeLSNumbers
    Interpolation(:,iCount) = pchip(configurationSetup.iLoudspeakerFreqFilter(iCount).freqVector,...
        configurationSetup.iLoudspeakerFreqFilter(iCount).freq,signal_to_play.freqVector);
end

% Native struct for frequency filter (replaces itaAudio freq-domain constructor)
frequencyFilter.freq = Interpolation;
frequencyFilter.samplingRate = iFs;

%% 
new_page = zeros(length(signal_to_play.time),length(activeLSNumbers));
scaler = zeros(1,length(activeLSNumbers));
stimulus = zeros(length(signal_to_play.time),length(activeLSNumbers));
signal_run = zeros(length(signal_to_play.time),length(activeLSNumbers));
levels = zeros(1,max(activeLSNumbers));

%% VBAP balance receive the correct SPL from the two speakers

angle_space = 360/length(configurationSetup.ls_dir); %
%Clockwise adjust to match odeon

if iAngles > 90 && iAngles <= 180 %Desired presentation to sequential LS
    if  iAngles <= 135
        %     s1 = 1; s2 = 19;
        s1 = activeLSNumbers(configurationSetup.ls_dir(:,1)==180);
        s2 = activeLSNumbers(configurationSetup.ls_dir(:,1)==90);
    else
        s1 = activeLSNumbers(configurationSetup.ls_dir(:,1)==90);
        s2 = activeLSNumbers(configurationSetup.ls_dir(:,1)==180);
    end
elseif iAngles > 180 && iAngles <= 270
    if  iAngles <= 225
        s1 = activeLSNumbers(configurationSetup.ls_dir(:,1)==180);
        s2 = activeLSNumbers(configurationSetup.ls_dir(:,1)==270);
    else
        s1 = activeLSNumbers(configurationSetup.ls_dir(:,1)==180);
        s2 = activeLSNumbers(configurationSetup.ls_dir(:,1)==270);
    end
elseif iAngles > 270 && iAngles <= 360
    if  iAngles <= 315
        s1 = activeLSNumbers(configurationSetup.ls_dir(:,1)==270);
        s2 = activeLSNumbers(configurationSetup.ls_dir(:,1)==0);
    else
    	s1 = activeLSNumbers(configurationSetup.ls_dir(:,1)==0);
        s2 = activeLSNumbers(configurationSetup.ls_dir(:,1)==270);
    end
elseif (iAngles >= 0 && iAngles <= 90)
    if  iAngles <= 45
        s1 = activeLSNumbers(configurationSetup.ls_dir(:,1)==0);
        s2 = activeLSNumbers(configurationSetup.ls_dir(:,1)==90);
    else
        s1 = activeLSNumbers(configurationSetup.ls_dir(:,1)==90);
        s2 = activeLSNumbers(configurationSetup.ls_dir(:,1)==0);
    end
end

%determine panning levels
%Coherent sum (6 dB) , so ^2 to calculate the correct energy level factor %Otherwise just sin/cos to pan law    
s1_level = (sin(abs(rem(iAngles,angle_space)/angle_space*(pi/2))))^2;
s2_level = (cos(abs(rem(iAngles,angle_space)/angle_space*(pi/2))))^2;

[s1_level, s2_level] = deal(max([s1_level s2_level]),min([s1_level  s2_level]));
if level ~= 'n'
    for idx = 1:length(level)
%         [s1(idx), s2(idx) s1_level(idx) s2_level(idx)] = equal_power_pan(iAngles(idx));
        if isinf(20*log10(10^((level(idx)/20))*s1_level(idx)))==1
            levels(s1(idx)) = 0;
        else
            levels(s1(idx)) = 20*log10(10^((level(idx)/20))*s1_level(idx));
        end
        if isinf(20*log10(10^((level(idx)/20))*s2_level(idx)))
            levels(s2(idx)) = 0;
        else
            levels(s2(idx)) = 20*log10(10^((level(idx)/20))*s2_level(idx));
        end
    end
end
% Native struct to hold filtered/leveled output (replaces itaAudio time-domain constructor)
signal_run_FILTER.time = zeros(length(signal_to_play.time), length(activeLSNumbers));
signal_run_FILTER.samplingRate = iFs;
signal_run_FILTER_LEVEL.time = zeros(length(signal_to_play.time), length(activeLSNumbers));
signal_run_FILTER_LEVEL.samplingRate = iFs;

%% What is happening here?
for idx = 1:length(activeLSNumbers)
    if isnan(signal_to_play.ch(idx).time)==1
       new_page(:,idx) = 0;
       fprintf('NaN Channel, take a look ch: %i\n',idx);
    else
      
        new_page(:,idx) = signal_to_play.ch(idx).time;
    end
    scaler(idx) = (sqrt(mean(new_page(:,idx).^2)))*2;
    if (scaler(idx)~=0) %Avoid NaN
        stimulus(:,idx) = new_page(:,idx)./scaler(idx);
    end
%     if level ~= 'n'
        signal_run(:,idx) = stimulus(:,idx).*repmat(10.^((levels(activeLSNumbers(idx)) - lsdBperVolt )./20),length(new_page(:,idx)),1);
%     else
%         signal_run(:,idx)= signal_to_play.ch(idx).time; 
%     end
    % Native spectral multiplication (replaces ita_multiply_spk)
    sigRunChFFT = fft(signal_run(:, idx));
    filterResp = frequencyFilter.freq(:, activeLSNumbers(idx));
    if length(filterResp) ~= length(sigRunChFFT)
        filterResp = interp1(linspace(0,1,length(filterResp)), filterResp, linspace(0,1,length(sigRunChFFT))).';
    end
    filtered = real(ifft(sigRunChFFT .* filterResp));
    signal_run_FILTER.time(:,idx) = filtered;
    signal_run_FILTER_LEVEL.time(:,idx) = signal_run_FILTER.time(:,idx) .* (Level_Factor(activeLSNumbers(idx)));
end
%%
calibrated_vbap = signal_run_FILTER_LEVEL;
% calibrated_vbap   = VBAP_DS*max(DSER.time); 
end

