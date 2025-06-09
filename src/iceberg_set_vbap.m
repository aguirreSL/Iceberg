function [signal_vbap] = iceberg_set_vbap(signal, DSER, iAngle, configurationSetup)
%% Define VBAP signal to the specific array using the most recent calibration
%
% [signalcalibrated] = iceberg_set_vbap(itaAudio,iAngle,configurationSetup)
% itaAudio each channel will be played at defined level and angle
% iAngle: Should be a vector as the same number of signal channels
% Author: Sergio
% Updated: 09/06/2025

%% Load values from the most recent calibration
iFs = signal.samplingRate;                          %Sample Frequency
signal = ita_convolve(signal, DSER);

%% Prepare Signal
        signal_vbap = itaAudio();        %Empty itaAudio object
        signal_vbap.samplingRate = iFs;  %set FS (if you use other than 44.1 this step is essential)
        signal_vbap.fftDegree = log2(iFs*signal.trackLength); 
        for nSignalsIndex = 1:signal.nChannels
            signal_vector(nSignalsIndex,:,:) = ring_VBAP(...
                iFs,signal.time(:,nSignalsIndex), iAngle(nSignalsIndex), configurationSetup);
            signal_ita(nSignalsIndex)  = itaAudio(squeeze(signal_vector(nSignalsIndex,:,:)), iFs, 'time');
            signal_vbap = ita_add(signal_vbap, signal_ita(nSignalsIndex));
        end

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
ls_groups = findLsPairs(configurationSetup.ls_dir(:,1)); 
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
