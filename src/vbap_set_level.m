function [signalcalibrated] = vbap_set_level(signal,level,iAngles,configurationSetup)
%% Define VBAP signal to the specific array using the most recent calibration
%
% [signalcalibrated] = vbap_set_level(itaAudio,level,iAngles,configurationSetup)
% itaAudio each channel will be played at defined level and angle
% level: Should be a vector as the same number of signal channels
% iAngles: Should be a vector as the same number of signal channels
% Author: Sergio
% Date: 21/11/21

%% Load values from the most recent calibration
Level_Factor = configurationSetup.new_Level_Factor; %LS Level adjust
iFs = signal.samplingRate;                          %Sample Frequency
lsdBperVolt = (20*log10((configurationSetup.iFactor)/2e-5)); %Factor
%% Prepare Signal
        signal_to_play = itaAudio();        %Empty itaAudio object
        signal_to_play.samplingRate = iFs;  %set FS (if you use other than 44.1 this step is essential)
        signal_to_play.fftDegree = log2(iFs*signal.trackLength); %At that point you got the FFTdegree to seconds already, right? t = (2^(fftDegree))/fs
        for nSignalsIndex = 1:signal.nChannels
            signal_vector(nSignalsIndex,:,:) = ring_VBAP(iFs,signal.time(:,nSignalsIndex),iAngles(nSignalsIndex),configurationSetup);
            signal_ita(nSignalsIndex)  = itaAudio(squeeze(signal_vector(nSignalsIndex,:,:)),iFs,'time');
            signal_to_play = ita_add(signal_to_play,signal_ita(nSignalsIndex));
        end

%% 
activeLSNumbers = zeros(1,length(configurationSetup.ls_dir));
for indexx= 1:length(configurationSetup.ls_dir)
iArrayP = find(configurationSetup.LSArray==configurationSetup.ls_dir(indexx,1),1);
activeLSNumbers(indexx) = iArrayP;
end
for iCount = activeLSNumbers
    Interpolation(:,iCount) = pchip(configurationSetup.iLoudspeakerFreqFilter(iCount).freqVector,...
        configurationSetup.iLoudspeakerFreqFilter(iCount).freq,signal_to_play.freqVector);
end

frequencyFilter = itaAudio(Interpolation,iFs,'freq');

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
signal_run_ita_FILTER       = itaAudio(zeros(length(signal_to_play.time),length(activeLSNumbers)),iFs,'time');
signal_run_ita_FILTER_LEVEL = itaAudio(zeros(length(signal_to_play.time),length(activeLSNumbers)),iFs,'time');

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
    signal_run_ita_dB = itaAudio(signal_run,iFs,'time');
    filtered = ita_multiply_spk(signal_run_ita_dB.ch(idx),frequencyFilter.ch(activeLSNumbers(idx)));
    signal_run_ita_FILTER.time(:,idx) = filtered.time;
    signal_run_ita_FILTER_LEVEL.time(:,idx) = signal_run_ita_FILTER.ch(idx).time.*(Level_Factor(activeLSNumbers(idx)));
end
%%
signalcalibrated = signal_run_ita_FILTER_LEVEL;


end

function [pansig,padsig,sig] = ring_VBAP(fs,iSignal,iAngle,configurationSetup) 
%% VBAP Archontis
% initial parameters (Here we can even move the sound source... I will
% update this to have as an option, now is static) .Sergio.
blocksize = fs/18.3750; % (~50msec)
if fs == 48000
    blocksize = fs/20; % (~50msec)
elseif fs == 96000
    blocksize = fs/40; % (~50msec)
end
hopsize = blocksize/2; % panning hopsize (update the panning values twice per bocksize) 
ls_num = length(configurationSetup.ls_dir);
sig = iSignal;
%% Parametros do sinal
Lsig = length(sig);
Nhop = ceil(Lsig/hopsize) + 2;
padsig = [zeros(hopsize,1); sig ;zeros(Nhop*hopsize - Lsig - hopsize,1)]; % zero padding
pansig = zeros(size(padsig,1), ls_num);
%% define the trajectory for panning (Quando sinal esta girando girando...)
static = ones(length((0:(Nhop-1)-1)'*(9*360)/(Nhop-1)),1); %verify here
iAngleCount = iAngle;
azis = iAngleCount*static;
eles = 0*azis;
% precompute VBAP matrix inversion (I LOVE this part)
ls_groups = findLsPairs(configurationSetup.ls_dir(:,1)); %it can return also a mesh for plotting
layoutInvMtx = invertLsMtx(configurationSetup.ls_dir(:,1), ls_groups);
%% do panning of signal to defined trajectory with overlap-add method
counter = 1;
window = hanning(blocksize);
spread = 0;

for idx = 0:hopsize:(Nhop-2)*hopsize
    winsig = padsig(idx+(1:blocksize),1).*window;
    azi = azis(counter);
    ele = eles(counter);
    gains   = vbap([azi ele], ls_groups, layoutInvMtx, spread); %Using VBAP
    gains2  = vbip([azi ele], ls_groups, layoutInvMtx, spread); %Using VBIP
    panwinsig = winsig*gains;
    pansig(idx+(1:blocksize),:) = pansig(idx+(1:blocksize),:) + panwinsig;
    counter = counter+1;
end
% truncate loudspeaker signals to original length (omit zeropadding)
pansig = pansig(hopsize+(1:Lsig),:);


end
%%
% 
% function [s1, s2, s1_level, s2_level] = equal_power_pan(angle)
% %EQUAL_POWER_PAN returns loudspeakers and scalars
% % Given a desired sound source angle, this function computes the two
% % nearest loudspeakers and the associated scalars to ensure the signal is
% % panned between them without a change in overall level
% %
% % Example:
% % >>[s1 s2 s1_level s2_level] = equal_power_pan(40)
% % s1 = 6
% % s2 = 5
% % s1_level = 0.8660
% % s2_level = 0.5000
% %
% % Then adjust your "signal" as follows:
% % s1_signal = signal.*s1_level
% % s2_signal = signal.*s2_level
% % Author: Owen
% % Date: 14/01/11
% 
% %establish loudspeaker locations and numbers
% % speaker_angles = 0:15:345;
% % speaker_nums = 1:25;
% speaker_angles = 0:90:345;
% speaker_nums = 1:6:24;
% 
% %fold over angles:
% angle = mod(angle+90,360);
% 
% %determine nearest loudspeaker to the desired angle
% [~,index] = sort(abs(speaker_angles - angle));
% 
% [s1, s2] = deal(speaker_nums(index(1)),speaker_nums(index(2)));
% 
% %fold over ring rotations:
% s1 = mod(s1-1,4)+1; s2 = mod(s2-1,4)+1;
% 
% %determine panning levels
% s1_level = sin(abs(rem(angle,15)/15*(pi/2)));
% s2_level = cos(abs(rem(angle,15)/15*(pi/2)));
% 
% %swap if s1 is quieter
% [s1_level, s2_level] = deal(max([s1_level s2_level]),min([s1_level  s2_level]));
% end
