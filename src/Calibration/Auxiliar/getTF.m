function [iLoudspeakerFreqFilter,rir,rirAdjusted] = getTF(iAverages,iFs,ifftDegree,inputChannel,iPlot)
%% Room Impulse Response Eriksholm Anechoic Chamber
% Microphone: B&K 1/2'' 4192 (Pressure field compensated) Vertically
% oriented, if you use a free field mic, orient towards each sound source
iLoudspeakerFreqFilter = itaResult();
% Measurement parameters
% iAverages = 2;
% ifftDegree = 15;
% iFs = 44100;
iStopMargin = 0.1;
iFreqRange = [20 24000];
iMode = 'rapid';
%Creating measurement
%Generate sweep 
%Frequency range 50 to 20,000 Hz
%Sample Frequency 44,100 Hz

%FFTDegree 16 | trackLength   = 1.486 s (t= (2^(fftDegree))/fs)
%FFTDegree 17 | trackLength   = 2.972 s (t= (2^(fftDegree))/fs)
%FFTDegree 18 | trackLength   = 5.944 s (t= (2^(fftDegree))/fs)

signal = ita_generate('swenlinsweep',[50 22627],0.1,iFs,ifftDegree);



%Measuring Loop
LS =1:24;


for iLoudspeaker = LS
    fprintf('Measuring Loudspeaker %i\n',iLoudspeaker)
    m = itaMSTF('freqRange', iFreqRange,...
        'pause', 1.5,...
        'inputChannels', inputChannel,...
        'outputChannels', 1,...
        'averages', iAverages,...
        'comment',['LS = ' num2str(iLoudspeaker)],...
        'outputamplification', -22,...
        'fftDegree',ifftDegree,...
        'stopMargin',iStopMargin,...
        'outputEqualizationFilters',[],...
        'excitation',signal);
    
    %Select Loudspeaker
    m.outputChannels = iLoudspeaker;
    %Select Signal
    m.excitation = signal;
    %Run measurement to the selected loudspeaker
    if strcmp(iMode,'complete')
        rirRaw(iLoudspeaker) = m.run_raw;
        rirLatency(iLoudspeaker) = m.run_latency;
        rirBackgroundNoise = m.run_backgroundNoise;
        rirSNR(iLoudspeaker) = m.run_snr;
        rir(iLoudspeaker) = m.run;
    else
        rir(iLoudspeaker) = m.run;
    end
    
    % Get the Room Impulse Response (RIR) in 
    % 1/3 octave|Filter|Smooth Freq|Normalize|Invert
    
    
    RIR_third_octave_inverted = ita_spk2frequencybands(...
        ita_invert_spk(...
        ita_normalize_spk(...
        ita_smooth_frequency(...
        ita_filter_bandpass(...
        rir(iLoudspeaker),...
        'upper',22627,'lower',20)...
        ))),'bandsperoctave',3,'freqRange',[50 22627]);
    
    RIR_third_octave_inverted.freqData(end) = RIR_third_octave_inverted.freqData(end-1);
    
    
    %Interpolation to smooth the filter
    Interpolation = pchip(RIR_third_octave_inverted.freqVector,RIR_third_octave_inverted.freq,signal.freqVector);
    %Transform filter in itaAudio
    magicFilter = (itaAudio(Interpolation,iFs,'freq'));
    %Aplying filter to the created measurement
    m.outputEqualizationFilters = magicFilter;
    %Save Filter coefs
    iLoudspeakerFreqFilter(iLoudspeaker) = RIR_third_octave_inverted;
    %Adjust Output
    m.outputamplification = -22;
    m.averages = 1;
    pause(2)
    %Run measurement again with adjusted individualized filter
    rirAdjusted(iLoudspeaker) = m.run;

    %Reset the measurement
    
    
    
end


%%
if iPlot == 1
    for iLoudspeaker = LS
        eval(['rir' num2str(iLoudspeaker) ' = ita_normalize_spk(ita_smooth_frequency(ita_filter_bandpass(rir(' num2str(iLoudspeaker) '),''lower'' , 20,''upper'', 22050)));']);
    end
    
    ita_plot_freq(merge(rir1,rir2,rir3,rir4,rir5,rir6,rir7,rir8,rir9,rir10,...
        rir11,rir12,rir13,rir14,rir15,rir16,rir17,rir18,rir19,rir20,rir21,...
        rir22,rir23,rir24))
    
    %%
    ylim([-40 20])
    xlim([45 17e3])
    set(gca,'fontsize',22)
    title(sprintf('Normalized Room Impulse Response.'))
    ylabel('Normalized Magnitude [dB]')
    xlabel('Frequency [Hz]')
    % Plot the boundaries of ± 3dB ..
    xu = [50 250 2000 16000]; %lower limit
    yu = [-7.67 -5 -5 -7.5];
    xo = [50 16000]; %upper limit
    yo = [3 3];
    plot(xu,yu,'Color',[0.4 0.4 0.4],'LineWidth',3,'LineStyle','--'); %draw boundaries
    plot(xo,yo,'Color',[0.4 0.4 0.4],'LineWidth',3,'LineStyle','--');
    % 16 kHz line and tolerance boundaries
    line([16e3 16e3], [-20 6],'LineWidth',3,'Color',[1 0 0],'LineStyle','--');
    line([50 50], [-20 6],'LineWidth',3,'Color',[1 0 0],'LineStyle','--');
    %%
    for iLoudspeaker = LS
        eval(['rir' num2str(iLoudspeaker) ' = ita_normalize_spk(ita_smooth_frequency(ita_filter_bandpass(rirAdjusted(' num2str(iLoudspeaker) '),''lower'' , 20,''upper'', 22050)));']);
    end
    
    ita_plot_freq(merge(rir1,rir2,rir3,rir4,rir5,rir6,rir7,rir8,rir9,rir10,...
        rir11,rir12,rir13,rir14,rir15,rir16,rir17,rir18,rir19,rir20,rir21,...
        rir22,rir23,rir24))
    
    
    ylim([-40 20])
    xlim([45 17e3])
    set(gca,'fontsize',22)
    title(sprintf('Normalized Room Impulse Response. (Adjusted)'))
    ylabel('Normalized Magnitude [dB]')
    xlabel('Frequency [Hz]')
    % Plot the boundaries of ± 3dB ..
    xu = [50 250 2000 16000]; %lower limit
    yu = [-7.67 -5 -5 -7.5];
    xo = [50 16000]; %upper limit
    yo = [3 3];
    plot(xu,yu,'Color',[0.4 0.4 0.4],'LineWidth',3,'LineStyle','--'); %draw boundaries
    plot(xo,yo,'Color',[0.4 0.4 0.4],'LineWidth',3,'LineStyle','--');
    % 16 kHz line and tolerance boundaries
    line([16e3 16e3], [-20 6],'LineWidth',3,'Color',[1 0 0],'LineStyle','--');
    line([50 50], [-20 6],'LineWidth',3,'Color',[1 0 0],'LineStyle','--');
end

end

