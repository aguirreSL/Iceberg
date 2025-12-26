%% Recording Backgroung Noise at Eriksholm Anechoic Room
% This script was used to measure the backgroung noise inside anechoic room
% Date created November 3rd, 2018. 
% Date last used September 3rd, 2020. 
% Settings:
% Microphone: B&K 1" 
% Preamp: B&K....
% Preamp source
% Sound card MOTU PCIe ...
% Acoustical Calibrator Brüel & Kjær 4231
%% Generate object to record via ITA-TB
recordObject = itaMSRecord();
recordObject.samplingRate=44100;
recordObject.inputChannels =1;
recordObject.freqRange = [11  11300];
recordObject.fftDegree=21;
recordObject.averages = 3;
recordObject.applyBandpass = 1;
%% Recording Calibrator Sound to reference 
% Remember to check the batteries
recordObject.comment = 'Calibrator 4231';%check calibrator model
signalFromCalibrator = recordObject.run;
% Calculate signal's rms from time domain
rmsSignal = rms(signalFromCalibrator.time); %signalFromCalibrator.rms (also works)
% evaluate in frequency domain at 1 kHz
valueAt1kHz = abs(signalFromCalibrator.freq2value(999.9481));
% Present in dB
valueAt1kHz_dB = 20*log10(abs(signalFromCalibrator.freq2value(999.9481)));
sprintf('  %.4f [VFS]\n%.4f [dBFS]',valueAt1kHz,valueAt1kHz_dB)

%% Check time domain
signalFromCalibrator.plot_time
hold on
plot(signalFromCalibrator.timeVector,ones(length(signalFromCalibrator.time),1)*rmsSignal,'LineWidth',4)
legend('Signal from calibrator [VFS]',sprintf('RMS value %.4f [VFS]',rmsSignal))
ylabel('Amplitude [VFS]')
%%
% Check frequency domain
signalFromCalibrator.plot_freq

disp('Remove calibrator, close the door and press any key to continue')
pause()

%% Factor 
vfs2splFactor = 1/valueAt1kHz; %1 [Pa] / valueAt1kHz [VFS] 

%% Checking 
adjustedSignal = signalFromCalibrator*vfs2splFactor;

%% SPL Reference 20 micro Pascals
referenceSPL = 2e-5;
adjustedSignalSPL = adjustedSignal/referenceSPL;
ita_plot_freq(adjustedSignalSPL)
ylabel('Sound Pressure Level [dB]');
title('Check Calibrator Level')
%% Recording the background noise (remember on removing the calibrator first!)
recordObject.comment = 'Lights on, Loudspeakers on';
signalFromAnechoicRoom = recordObject.run;
%% Plot All on linear
filterAllOn = ita_filter_bandpass(signalFromAnechoicRoom,'lower',11,'upper',11313);
filterAllOn_adjusted = ita_spk2level(filterAllOn*vfs2splFactor,1,'added');
filterAllOn_adjusted_averaged = ita_spk2level((filterAllOn*vfs2splFactor)/referenceSPL,0,'added');
filterAllOn_adjusted_SPL = filterAllOn_adjusted/referenceSPL; 
auxVector = filterAllOn_adjusted_SPL.freqData_dB;
figure();
set(gcf,'units','normalized','outerposition',[0 0 1 1])
bar(1:10,auxVector(1:10))
set(gca,'XTick',1:10,'XTicklabel',[16 31.5 62.5 125 250 500 1000 2000 4000 8000])
set(gca,'yscale','lin','YMinorTick','off')
ylim([-30 50])
title(recordObject.comment )
xlim([0.5 10.5])
ylabel('Sound Pressure Level [dB]')
xlabel('Frequency [Hz]')
title(sprintf('SPL_{eq} %.2f dB | Lights and Loudspeakers on ',filterAllOn_adjusted_averaged.freqData_dB))
set(gca,'FontSize',18)

%% Plot All on A weight
signalFromAnechoicRoom_ON_A = ita_mpb_filter(signalFromAnechoicRoom,'a-weight');
filterAllOn_A = ita_filter_bandpass(signalFromAnechoicRoom_ON_A,'lower',11,'upper',11313);
filterAllOn_adjusted_A = ita_spk2level(filterAllOn_A*vfs2splFactor,1,'added');
filterAllOn_adjusted_averaged_A = ita_spk2level((filterAllOn_A*vfs2splFactor)/referenceSPL,0,'added');
filterAllOn_adjusted_SPL_A = filterAllOn_adjusted_A/referenceSPL; 
auxVector = filterAllOn_adjusted_SPL_A.freqData_dB;
figure();
set(gcf,'units','normalized','outerposition',[0 0 1 1])
bar(1:10,auxVector(1:10))
set(gca,'XTick',1:10,'XTicklabel',[16 31.5 62.5 125 250 500 1000 2000 4000 8000])
set(gca,'yscale','lin','YMinorTick','off')
ylim([-30 50])
title(recordObject.comment )
xlim([0.5 10.5])
ylabel('Sound Pressure Level [dBA]')
xlabel('Frequency [Hz]')
title(sprintf('SPL(A)_{eq} %.2f dB | Lights and Loudspeakers on ',filterAllOn_adjusted_averaged_A.freqData_dB))
set(gca,'FontSize',18)

%% All off
recordObject.comment = 'Lights off, Loudspeakers off';
signalFromAnechoicRoom_OFF = recordObject.run;
%% Plot it
filterAllOff = ita_filter_bandpass(signalFromAnechoicRoom_OFF,'lower',11,'upper',11313);
Filtersignal_off_adjusted = ita_spk2level(filterAllOff*vfs2splFactor,1,'added');
Filtersignal_off_adjusted_averaged = ita_spk2level((filterAllOff*vfs2splFactor)/referenceSPL,0,'added');
filterAllOff_adjusted_SPL = Filtersignal_off_adjusted/referenceSPL;
vec = filterAllOff_adjusted_SPL.freqData_dB;
figure();
set(gcf,'units','normalized','outerposition',[0 0 1 1])
bar(1:10,vec(1:10))
set(gca,'XTick',1:10,'XTicklabel',[16 31.5 62.5 125 250 500 1000 2000 4000 8000])
set(gca,'yscale','lin','YMinorTick','off')
ylim([-30 50])
title(recordObject.comment )
xlim([0.5 10.5])
ylabel('Sound Pressure Level [dB]')
xlabel('Frequency [Hz]')
title(sprintf('SPL_{eq} %.2f dB | Lights and Loudspeakers off ',Filtersignal_off_adjusted_averaged.freqData_dB))
set(gca,'FontSize',18)

%% ALL off A wheight
Filter_A_off = ita_mpb_filter(signalFromAnechoicRoom_OFF,'a-weight')
Filtersignal_off_A = ita_filter_bandpass(Filter_A_off,'lower',11,'upper',11313)
Filtersignal_off_A_adjusted = ita_spk2level(Filtersignal_off_A*vfs2splFactor,1,'added');
Filtersignal_off_SPL_A_adjusted = Filtersignal_off_A_adjusted/referenceSPL;
Filtersignal_off_A_adjusted_averaged = ita_spk2level((Filtersignal_off_A*vfs2splFactor)/referenceSPL,0,'added');
Filtersignal_off_A_adjusted_averaged.freqData_dB;
vec = Filtersignal_off_SPL_A_adjusted.freqData_dB;
figure();
set(gcf,'units','normalized','outerposition',[0 0 1 1])
bar(1:10,vec(1:10))
set(gca,'XTick',1:10,'XTicklabel',[16 31.5 62.5 125 250 500 1000 2000 4000 8000])
set(gca,'yscale','lin','YMinorTick','off')
ylim([-30 40])
title(sprintf('SPL(A)_{eq} %.2f dBA | Lights and Loudspeakers off ',Filtersignal_off_A_adjusted_averaged.freqData_dB))
xlim([0.5 10.5])
ylabel('Sound Pressure Level [dBA]')
xlabel('Frequency [Hz]')
set(gca,'FontSize',18)
%%

