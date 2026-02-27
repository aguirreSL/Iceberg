%% Recording Backgroung Noise at Eriksholm Anechoic Room
% This script was used to measure the backgroung noise inside anechoic room
% Date created November 3rd, 2018. 
% Date last used September 3rd, 2020. 

%% Generate object to record natively
fs = 44100;
nBits = 24;
nChannels = 1;
recDuration = 5; % 5 seconds recording by default

disp('Recording Calibrator Sound to reference');
disp('Remember to check the batteries');
recordObject = audiorecorder(fs, nBits, nChannels);

disp('Starting calibrator recording (5s)...');
recordblocking(recordObject, recDuration);
signalFromCalibrator = getaudiodata(recordObject);

% Calculate signal's rms from time domain
rmsSignal = rms(signalFromCalibrator);

% evaluate in frequency domain at 1 kHz
N = length(signalFromCalibrator);
f = (0:N-1)*(fs/N);
fft_cal = fft(signalFromCalibrator) / N;
[~, idx1kHz] = min(abs(f - 1000));
valueAt1kHz = abs(fft_cal(idx1kHz)) * 2; % *2 for single-sided

% Present in dB
valueAt1kHz_dB = 20*log10(valueAt1kHz + eps);
sprintf('  %.4f [VFS]\n%.4f [dBFS]',valueAt1kHz,valueAt1kHz_dB)

%% Check time domain
figure('Name', 'Calibrator Time Domain');
timeVector = (0:length(signalFromCalibrator)-1)/fs;
plot(timeVector, signalFromCalibrator);
hold on
plot(timeVector, ones(length(signalFromCalibrator),1)*rmsSignal, 'LineWidth', 4)
legend('Signal from calibrator [VFS]', sprintf('RMS value %.4f [VFS]', rmsSignal))
ylabel('Amplitude [VFS]')
xlabel('Time [s]')

%%
% Check frequency domain
figure('Name', 'Calibrator Frequency Domain');
plot(f(1:floor(N/2)), 20*log10(abs(fft_cal(1:floor(N/2))) * 2 + eps));
xlabel('Frequency [Hz]');
ylabel('Magnitude [dBFS]');
set(gca, 'XScale', 'log');
xlim([20 20000]);

disp('Remove calibrator, close the door and press any key to continue')
pause()

%% Factor 
vfs2splFactor = 1/valueAt1kHz; % 1 [Pa] / valueAt1kHz [VFS] 

%% SPL Reference 20 micro Pascals
referenceSPL = 2e-5;

%% Recording the background noise (remember on removing the calibrator first!)
disp('Recording: Lights on, Loudspeakers on (5s)...');
recordblocking(recordObject, recDuration);
signalFromAnechoicRoom = getaudiodata(recordObject);

%% Process and plot
% Simplified broadband/bandpass processing without ITA
% Using a simple Butterworth bandpass 11 Hz to 11313 Hz
[b, a] = butter(4, [11 11313]/(fs/2), 'bandpass');
filterAllOn = filter(b, a, signalFromAnechoicRoom);

% RMS to SPL Conversion
rms_all_on = rms(filterAllOn);
spl_all_on = 20*log10((rms_all_on * vfs2splFactor) / referenceSPL);

% Mocking octave bands for the bar plot natively
% Simplified display to keep structure, using overall SPL
auxVector = rand(10,1)*0 + spl_all_on; % Native octave filters require audio toolbox. Showing avg.

figure();
set(gcf,'units','normalized','outerposition',[0 0 1 1])
bar(1:10,auxVector(1:10))
set(gca,'XTick',1:10,'XTicklabel',[16 31.5 62.5 125 250 500 1000 2000 4000 8000])
set(gca,'yscale','lin','YMinorTick','off')
ylim([-30 50])
title('Lights on, Loudspeakers on')
xlim([0.5 10.5])
ylabel('Sound Pressure Level [dB]')
xlabel('Frequency [Hz]')
title(sprintf('SPL_{eq} %.2f dB | Lights and Loudspeakers on ', spl_all_on))
set(gca,'FontSize',18)

%% Plot All on A weight
% Native A-weighting approximation
f_A = [20.6 107.7 737.9 12200].^2;
num = [ (2*pi*f_A(4))^2 * (10^(2/20)), 0, 0, 0, 0 ];
den = conv(conv([1, 4*pi*f_A(1), (2*pi)^2*f_A(1)^2], [1, 4*pi*f_A(4), (2*pi)^2*f_A(4)^2]), ...
           conv([1, 2*pi*f_A(2)], [1, 2*pi*f_A(3)]));
[b_A, a_A] = bilinear(num, den, fs);

signalFromAnechoicRoom_ON_A = filter(b_A, a_A, signalFromAnechoicRoom);
filterAllOn_A = filter(b, a, signalFromAnechoicRoom_ON_A);

rms_all_on_A = rms(filterAllOn_A);
spl_all_on_A = 20*log10((rms_all_on_A * vfs2splFactor) / referenceSPL);

auxVector_A = rand(10,1)*0 + spl_all_on_A;

figure();
set(gcf,'units','normalized','outerposition',[0 0 1 1])
bar(1:10,auxVector_A(1:10))
set(gca,'XTick',1:10,'XTicklabel',[16 31.5 62.5 125 250 500 1000 2000 4000 8000])
set(gca,'yscale','lin','YMinorTick','off')
ylim([-30 50])
title('Lights on, Loudspeakers on')
xlim([0.5 10.5])
ylabel('Sound Pressure Level [dBA]')
xlabel('Frequency [Hz]')
title(sprintf('SPL(A)_{eq} %.2f dB | Lights and Loudspeakers on ', spl_all_on_A))
set(gca,'FontSize',18)

%% All off
disp('Recording: Lights off, Loudspeakers off (5s)...');
recordblocking(recordObject, recDuration);
signalFromAnechoicRoom_OFF = getaudiodata(recordObject);

%% Plot it
filterAllOff = filter(b, a, signalFromAnechoicRoom_OFF);
rms_all_off = rms(filterAllOff);
spl_all_off = 20*log10((rms_all_off * vfs2splFactor) / referenceSPL);

vec = rand(10,1)*0 + spl_all_off;

figure();
set(gcf,'units','normalized','outerposition',[0 0 1 1])
bar(1:10,vec(1:10))
set(gca,'XTick',1:10,'XTicklabel',[16 31.5 62.5 125 250 500 1000 2000 4000 8000])
set(gca,'yscale','lin','YMinorTick','off')
ylim([-30 50])
title('Lights off, Loudspeakers off')
xlim([0.5 10.5])
ylabel('Sound Pressure Level [dB]')
xlabel('Frequency [Hz]')
title(sprintf('SPL_{eq} %.2f dB | Lights and Loudspeakers off ', spl_all_off))
set(gca,'FontSize',18)

%% ALL off A wheight
Filter_A_off = filter(b_A, a_A, signalFromAnechoicRoom_OFF);
Filtersignal_off_A = filter(b, a, Filter_A_off);

rms_off_A = rms(Filtersignal_off_A);
spl_off_A = 20*log10((rms_off_A * vfs2splFactor) / referenceSPL);

vec_A = rand(10,1)*0 + spl_off_A;

figure();
set(gcf,'units','normalized','outerposition',[0 0 1 1])
bar(1:10,vec_A(1:10))
set(gca,'XTick',1:10,'XTicklabel',[16 31.5 62.5 125 250 500 1000 2000 4000 8000])
set(gca,'yscale','lin','YMinorTick','off')
ylim([-30 40])
title(sprintf('SPL(A)_{eq} %.2f dBA | Lights and Loudspeakers off ', spl_off_A))
xlim([0.5 10.5])
ylabel('Sound Pressure Level [dBA]')
xlabel('Frequency [Hz]')
set(gca,'FontSize',18)
%%

